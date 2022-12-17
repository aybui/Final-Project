import java.io.*;
import java.io.File;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Arrays;

public class TranslateHomologs3 {

	public static void main(String[] args) {

		String inFileStem = args[0]; 
		File refFile = new File( args[1] ); 
		int nLoci = Integer.parseInt(args[2]); 
		int K = 15;

		ArrayList<String> refNames = readReferenceFile( refFile );
		for( int locus = 1; locus <= nLoci; locus++) {

			System.out.print("\rL"+locus);

			// read in homolog file for this locus
			ArrayList<Sequence> sequences = readHomologFile(new File(inFileStem+locus+".fasta"));

			// put reference seqs in a list
			ArrayList<Sequence> references = new ArrayList<>();
			for(Sequence sequence : sequences){
				for(String refName : refNames){
					if(sequence.getHeader().startsWith(">"+refName)){
						references.add(sequence);
						break;
					}
				}
			}

			// put homologs in a list
			ArrayList<Sequence> homologs = new ArrayList<>();
			for(Sequence sequence : sequences) {
				if( ! references.contains(sequence) ) {
					homologs.add(sequence);
				}
			}

			ArrayList<TranslatedSequence> translatedHomologs = new ArrayList<>();

			for(Sequence homolog : homologs)
			{
				for(Sequence reference : references) {

					String refSeq = getBestReadingFrame(reference.getSeq());
					if(refSeq == null) {
						continue;
					}

					TranslatedSequence translatedHomolog = translateHomologInAllReadingFrames(homolog, hashKmers(refSeq,K), K);

					if(translatedHomolog.wasSuccessful()) {
						translatedHomologs.add(translatedHomolog);
						break;
					}
				}

			}

			String outFileStem = inFileStem.split("/")[inFileStem.split("/").length-1];
			writeFiles(translatedHomologs, outFileStem, locus);
		}

		System.out.println();
		System.out.println("Done!");

	}

	static void writeFiles(ArrayList<TranslatedSequence> list, String outFileStem, int locus) {
		try {
			new File("../HomologsAA/").mkdir();
			new File("../HomologsNT12/").mkdir();
			new File("../HomologsNT123/").mkdir();
			new File("../HomologsNTALL/").mkdir();

			BufferedWriter bwAA = new BufferedWriter( new FileWriter( new File( "../HomologsAA/" + outFileStem + locus + ".fasta")));
			BufferedWriter bwNT12 = new BufferedWriter( new FileWriter( new File( "../HomologsNT12/" + outFileStem + locus + ".fasta")));
			BufferedWriter bwNT123 = new BufferedWriter( new FileWriter( new File( "../HomologsNT123/" + outFileStem + locus + ".fasta")));
			BufferedWriter bwNTALL = new BufferedWriter( new FileWriter( new File( "../HomologsNTALL/" + outFileStem + locus + ".fasta")));

			for(TranslatedSequence tseq : list) {
				bwAA.write(tseq.getHomologAA().getHeader() + "\n" + tseq.getHomologAA().getSeq() + "\n");
				bwNT12.write(tseq.getHomologNT12().getHeader() + "\n" + tseq.getHomologNT12().getSeq() + "\n");
				bwNT123.write(tseq.getHomologNT123().getHeader() + "\n" + tseq.getHomologNT123().getSeq() + "\n");
				bwNTALL.write(tseq.getHomologNTALL().getHeader() + "\n" + tseq.getHomologNTALL().getSeq() + "\n");
			}

			bwAA.flush();
			bwNT12.flush();
			bwNT123.flush();
			bwNTALL.flush();

			bwAA.close();
			bwNT12.close();
			bwNT123.close();
			bwNTALL.close();
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}

	static TranslatedSequence translateHomologInAllReadingFrames(Sequence homolog, HashMap<String,Integer> map, int K) {
		TranslatedSequence tseq = new TranslatedSequence();

		int[] refPos = new int[homolog.getSeq().length()];
		int bestNAAs = 0;
		String tempAA = "";

		for(int f = 0; f < 6; f++) {
			if(f<3){	//forward reading frames
				tempAA = translate(homolog.getSeq().substring(f));
			}else{		//reverse reading frames
				tempAA = translate(getRevComp(homolog.getSeq()).substring(f-3));
			}

			int lenAA = tempAA.length();
			Arrays.fill(refPos, -1);

			//position in translated amino acid sequence
			for(int j=0; j<lenAA-K+1; j++){
				String kmer=tempAA.substring(j,j+K);
				Integer value = map.get(kmer);
				if(value!=null){
					for(int x=0; x<K; x++){
						refPos[j+x]=value+x;
					}
				}
			}

			//fill in obvious ref positions that are missing because of a mismatch
			char posType[] = new char[lenAA];
			int refPos2[] = new int[lenAA];

			for(int x=0; x<lenAA; x++){
				if(refPos[x]==-1){
					posType[x]='-';
				}else{
					posType[x]='+';
				}
				if(x>0 && x<lenAA-1){
					if(refPos[x-1]==-1 && refPos[x]>-1){  //
						posType[x]='B';
					}else if(refPos[x]>-1 && refPos[x+1]==-1){
						posType[x]='E';
					}else{

					}
				}
			}

			for(int x=0; x<lenAA; x++){
				if(posType[x]=='B'){
					//find previous beginning and check difference
					for(int z=x-1; z>=0; z--){
						if(posType[z]=='E'){
							if(z-x == refPos[z]-refPos[x]){
								for(int w=z; w<=x; w++){
									refPos[w]=refPos[z]+(w-z);
								}
							}
							break;
						}
					}
				}
			}

			//count the number of valid amino acids that were translated for this reading frame
			int currNAAs=0;
			for(int x=0; x<lenAA; x++){
				if(refPos[x]>-1){currNAAs++;}
			}

			if(currNAAs>bestNAAs){  //found better reading frame...store sequences
				bestNAAs=currNAAs;
				//isolate the valid amino acid sequence and corresponding nucleotide sequences
				String homologAA = "";
				String homologNT123 = "";
				String homologNT12 = "";
				String homologNTALL = homolog.getSeq();
				for(int x=0; x<lenAA; x++){
					if(refPos[x]>-1){
						homologAA+=tempAA.charAt(x);
						homologNT12+=homolog.getSeq().substring((f%3)+x*3,(f%3)+x*3+2);
						homologNT123+=homolog.getSeq().substring((f%3)+x*3,(f%3)+x*3+3);
					}
				}

				tseq.setHomologAA(new Sequence(homolog.header, homologAA));
				tseq.setHomologNT12(new Sequence(homolog.header, homologNT12));
				tseq.setHomologNT123(new Sequence(homolog.header, homologNT123));
				tseq.setHomologNTALL(new Sequence(homolog.header, homologNTALL));
			}
		}
		if(bestNAAs > 0) { //success
			tseq.wasSuccessful(true);
		}

		return tseq;
	}

	static HashMap<String,Integer> hashKmers(String seq, int K){
		HashMap<String,Integer> map = new HashMap<>();

		for(int i = 0; i < seq.length()-K+1; i++) {
			String kmer = seq.substring(i,i+K);

			Integer value = map.get(kmer);
			if(value==null) {
				map.put(kmer, (Integer)i); //value is position in reference amino acid sequence at which the kmer starts
			}
			else {
				map.put(kmer, (Integer)(-9999));
			}
		}

		return map;
	}

	static ArrayList<Sequence> readHomologFile(File homologFile){

		ArrayList<Sequence> list = new ArrayList<>();

		try{
			BufferedReader br = new BufferedReader(new FileReader(homologFile));
			String header = null;
			String seq = null;

			while( (header = br.readLine()) != null && (seq = br.readLine()) != null) {
				if(seq.length() >= 3) {
					list.add(new Sequence(header, seq));
				}
			}

		} catch(IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}

		return list;
	}

	static String getBestReadingFrame(String seq) {

		int numStops = 0;
		String temp = null;
		for(int i = 0; i < 6; i++) {

			if( i < 3 ){ //forward reading frames
				temp = translate(seq.substring(i));
			}
			else{ //reverse reading frames
				temp = translate(getRevComp(seq).substring(i-3));
			}

			numStops = 0;

			//count stop codons
			for (char c : temp.toCharArray()){
				if(c == '*'){
					numStops++;
				}
			}
			if(numStops <= 1){
				break;
			}
		}
		if(numStops > 1){
			temp = null;
		}

		return temp;
	}

	static ArrayList<String> readReferenceFile(File refFile) {

		ArrayList<String> list = new ArrayList<>();
		String line = null;

		try{
			BufferedReader br = new BufferedReader( new FileReader( refFile ) );

			while( (line = br.readLine()) != null)
			{
				list.add(line);
			}
		}
		catch(IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}

		return list;
	}

	static String translate(String ntSeq){
		String aaSeq="";
		char b1,b2,b3;
		for(int i=0; i<ntSeq.length()-2; i+=3){
			b1=ntSeq.charAt(i);
			b2=ntSeq.charAt(i+1);
			b3=ntSeq.charAt(i+2);

			switch(b1){
				case 'A':
					switch(b2){
						case 'A':
							switch(b3){
								case 'A': aaSeq+='K'; break;
								case 'T': aaSeq+='N'; break;
								case 'C': aaSeq+='N'; break;
								case 'G': aaSeq+='K'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'T':
							switch(b3){
								case 'A': aaSeq+='I'; break;
								case 'T': aaSeq+='I'; break;
								case 'C': aaSeq+='I'; break;
								case 'G': aaSeq+='M'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'C':
							switch(b3){
								case 'A': aaSeq+='T'; break;
								case 'T': aaSeq+='T'; break;
								case 'C': aaSeq+='T'; break;
								case 'G': aaSeq+='T'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'G':
							switch(b3){
								case 'A': aaSeq+='R'; break;
								case 'T': aaSeq+='S'; break;
								case 'C': aaSeq+='S'; break;
								case 'G': aaSeq+='R'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						default:  aaSeq+='X'; break;
					}
					break;
				case 'T':
					switch(b2){
						case 'A':
							switch(b3){
								case 'A': aaSeq+='*'; break;
								case 'T': aaSeq+='Y'; break;
								case 'C': aaSeq+='Y'; break;
								case 'G': aaSeq+='*'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'T':
							switch(b3){
								case 'A': aaSeq+='L'; break;
								case 'T': aaSeq+='F'; break;
								case 'C': aaSeq+='F'; break;
								case 'G': aaSeq+='L'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'C':
							switch(b3){
								case 'A': aaSeq+='S'; break;
								case 'T': aaSeq+='S'; break;
								case 'C': aaSeq+='S'; break;
								case 'G': aaSeq+='S'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'G':
							switch(b3){
								case 'A': aaSeq+='*'; break;
								case 'T': aaSeq+='C'; break;
								case 'C': aaSeq+='C'; break;
								case 'G': aaSeq+='W'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						default:  aaSeq+='X'; break;
					}
					break;
				case 'C':
					switch(b2){
						case 'A':
							switch(b3){
								case 'A': aaSeq+='Q'; break;
								case 'T': aaSeq+='H'; break;
								case 'C': aaSeq+='H'; break;
								case 'G': aaSeq+='Q'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'T':
							switch(b3){
								case 'A': aaSeq+='L'; break;
								case 'T': aaSeq+='L'; break;
								case 'C': aaSeq+='L'; break;
								case 'G': aaSeq+='L'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'C':
							switch(b3){
								case 'A': aaSeq+='P'; break;
								case 'T': aaSeq+='P'; break;
								case 'C': aaSeq+='P'; break;
								case 'G': aaSeq+='P'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'G':
							switch(b3){
								case 'A': aaSeq+='R'; break;
								case 'T': aaSeq+='R'; break;
								case 'C': aaSeq+='R'; break;
								case 'G': aaSeq+='R'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						default:  aaSeq+='X'; break;
					}
					break;
				case 'G':
					switch(b2){
						case 'A':
							switch(b3){
								case 'A': aaSeq+='E'; break;
								case 'T': aaSeq+='D'; break;
								case 'C': aaSeq+='D'; break;
								case 'G': aaSeq+='E'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'T':
							switch(b3){
								case 'A': aaSeq+='V'; break;
								case 'T': aaSeq+='V'; break;
								case 'C': aaSeq+='V'; break;
								case 'G': aaSeq+='V'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						case 'C':
							switch(b3){
								case 'A': aaSeq+='A'; break;
								case 'T': aaSeq+='A'; break;
								case 'C': aaSeq+='A'; break;
								case 'G': aaSeq+='A'; break;
								default:  aaSeq+='A'; break;
							}
							break;
						case 'G':
							switch(b3){
								case 'A': aaSeq+='G'; break;
								case 'T': aaSeq+='G'; break;
								case 'C': aaSeq+='G'; break;
								case 'G': aaSeq+='G'; break;
								default:  aaSeq+='X'; break;
							}
							break;
						default:  aaSeq+='X'; break;
					}
					break;
				default:  aaSeq+='X'; break;
			}
		}
		return aaSeq;
	}

	static String getRevComp(String seq) {

		String revComp = "";
		int len = seq.length();
		for (int i = 0; i < len; i++) {
			revComp += seq.charAt(len - (i + 1));
		}
		if (seq.length() != revComp.length()) {
			System.out.println("ERROR 1 IN getRevCOMP!!!");
		}

		revComp = revComp.replace('A', 'Z');
		revComp = revComp.replace('T', 'A');
		revComp = revComp.replace('Z', 'T');

		revComp = revComp.replace('C', 'Q');
		revComp = revComp.replace('G', 'C');
		revComp = revComp.replace('Q', 'G');

		//REPLACE ALL OTHER BASES WITH 'N'
		for (int i = 0; i < len; i++) {
			char c = revComp.charAt(i);
			if (c != '-' && c != 'A' && c != 'T' && c != 'C' && c != 'G') {
				if (len > i + 1) {
					revComp = revComp.substring(0, i) + 'N' + revComp.substring(i + 1);
				} else {
					revComp = revComp.substring(0, i) + 'N';
				}
			}
		}
		if (revComp.length() != len) {
			System.out.println("ERROR 2 IN getRevCOMP!!!");
		}

		return revComp;
	}

	private static class Sequence {

		private String header;
		private String seq;

		public Sequence(String header, String seq){
			this.header = header;
			this.seq = seq;
		}

		public String getHeader() {
			return header;
		}

		public String getSeq() {
			return seq;
		}
	}

	private static class TranslatedSequence {
		private Sequence homologAA;
		private Sequence homologNT12;
		private Sequence homologNT123;
		private boolean successfulTranslation;
		Sequence HomologNTALL;

		public Sequence getHomologAA() {
			return homologAA;
		}

		public void setHomologAA(Sequence homologAA) {
			this.homologAA = homologAA;
		}

		public Sequence getHomologNT12() {
			return homologNT12;
		}

		public void setHomologNT12(Sequence homologNT12) {
			this.homologNT12 = homologNT12;
		}

		public Sequence getHomologNT123() {
			return homologNT123;
		}

		public void setHomologNT123(Sequence homologNT123) {
			this.homologNT123 = homologNT123;
		}

		public Sequence getHomologNTALL() {
			return HomologNTALL;
		}

		public void setHomologNTALL(Sequence homologNTALL) {
			HomologNTALL = homologNTALL;
		}

		public void wasSuccessful(boolean bool) {
			this.successfulTranslation = bool;
		}

		public boolean wasSuccessful() {
			return successfulTranslation;
		}
	}

}
