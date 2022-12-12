import java.io.*;
import java.util.Date;
import java.util.Arrays;
//ARL
//VERSION 3 CREATES SCRIPTS THAT DO NOT HAVE OUTGROUPS.
//Feb 22 2014  //THIS PROGRAM SETS UP THE FILE NEED TO GET THE RAxML ANALYSES GOING
//Mar 18 2015  //THIS UPDATE ADDS A FOLDERTAIL VARIABLE 
public class SetupRAxML3{

  public static void main(String[] args){
      try{

      	String project=args[0];
		int nLoci=Integer.parseInt(args[1]);			//maximum locus number
		int nThreadsA=Integer.parseInt(args[2]);		//threads for independent runs
		int nThreadsB=Integer.parseInt(args[3]);		//threads for concatenated run
		int nBootRepsA=Integer.parseInt(args[4]);		//boot reps for independent runs
		int nBootRepsB=Integer.parseInt(args[5]);		//boot reps for concatenated run			
		String folderTail="";
		if(args.length>6){folderTail=args[6];}			//e.g. _2alleles or _conSeqs or _mixturePloidy

		BufferedReader br=null; 
		
		//first determine the number of taxa ASSUME THAT FIRST file CONTAINS FULL NUMBER OF TAXA (EVEN IF SEQUENCES ARE EMPTY)
		int nTaxa=0;
		for(int loc=1; loc<=nLoci; loc++){
			if(!new File("../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta").exists()){ continue; }	
			
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   "../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta" ) ));
			String tempS=br.readLine();
			while(tempS!=null){
				if(tempS.startsWith(">")){
					nTaxa++;
				}
				tempS=br.readLine();
			}
			br.close();
			break;
		}	


		//now determine which sequences are missing
		boolean seqMissing[][] = new boolean[nTaxa][nLoci];
		int nSites[][] = new int[nTaxa][nLoci];
		for(int loc=1; loc<=nLoci; loc++){
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   "../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta" ) ));
			if(!new File("../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta").exists()){ continue; }	
			for(int i=0; i<nTaxa; i++){
				String tempS=br.readLine().replace(" ","_").replace("#","_").replace(".","_").replace("-","_");	//header
				tempS=br.readLine().toUpperCase();	//sequence
				nSites[i][loc-1]=tempS.length();
				if(tempS.replace("N","").replace("-","").length()==0){
					seqMissing[i][loc-1]=true;
				}
			}
			br.close();
		}
		br.close();
			
		//for each locus determine the number of taxa, and vice versa
		int nTaxaEach[] = new int[nLoci];
		int nLociEach[] = new int[nTaxa];
		int nSitesEach[] = new int[nLoci];
		
		for(int loc=1; loc<=nLoci; loc++){
			for(int i=0; i<nTaxa; i++){
				if(!seqMissing[i][loc-1]){
					nTaxaEach[loc-1]++;
					nLociEach[i]++;
					nSitesEach[loc-1]=Math.max(nSitesEach[loc-1],nSites[i][loc-1]);
				}
			}
		}
		
		//write the sequences in phylip format
		new File("../RAxML"+folderTail).mkdir();

		BufferedWriter bwXX = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "../RAxML"+folderTail+"/"+project+"_ConcatLoci.sh" )) ));
		bwXX.write("raxmlHPC -T "+nThreadsB+" -s "+project+"_ConcatLoci.phylip -q "+project+"_ConcatLoci.charSets -n "+project+"_ConcatLoci -m GTRGAMMA -p "+(int)(Math.random()*100000)+" -x "+(int)(Math.random()*100000)+" -# "+nBootRepsB+" -f a "+"\n");
		bwXX.flush();
		bwXX.close();

		BufferedWriter bwX = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "../RAxML"+folderTail+"/"+project+"_SeparateLoci.sh" )) ));

		BufferedWriter bwCharSets = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "../RAxML"+folderTail+"/"+project+"_ConcatLoci.charSets" )) ));

		BufferedWriter bwTEMP = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "temp" )) ));

		int lastCharPos=0;
		boolean firstGoodLocus=true;
		int nGoodLoci=0;
		int nGoodTaxa=0;
		
		for(int loc=1; loc<=nLoci; loc++){
			System.out.print("Locus "+loc);
			if(!new File("../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta").exists()){ System.out.println(" No File"); continue; }				

			if(nTaxaEach[loc-1]==0 || nSitesEach[loc-1]==0){System.out.println(" Empty matrix"); continue;}	//bail out if the matrix is empty

			System.out.println();

			bwCharSets.write("DNA,  L"+loc+"  =  "+(lastCharPos+1)+"-");
			lastCharPos+=nSitesEach[loc-1];
			bwCharSets.write((lastCharPos)+"\n");

			bwX.write("raxmlHPC -T "+nThreadsA+" -s "+project+"_L"+loc+".phylip -n "+project+"_L"+loc+" -m GTRGAMMA -p "+(int)(Math.random()*100000)+" -x "+(int)(Math.random()*100000)+" -# "+nBootRepsA+" -f a ");
			bwX.write("\n");

			nGoodLoci++;
			
			//this time write a phylip format alignment
			BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "../RAxML"+folderTail+"/"+project+"_L"+loc+".phylip" )) ));
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   "../TrimmedAlignments"+folderTail+"/"+project+"_L"+loc+".fasta" ) ));
			bw.write(nTaxaEach[loc-1]+" "+nSitesEach[loc-1]+"\n");
			String seqName;
			String seq;

			for(int i=0; i<nTaxa; i++){
				seqName=br.readLine();
				seq=br.readLine();
				seqName=seqName.replace(" ","_").replace("#","_").replace(".","_").replace("-","_");
				if(!seqMissing[i][loc-1]){
					bw.write(seqName.substring(1)+"\t"+seq+"\n");					
				}
				if(firstGoodLocus && nLociEach[i]>0){
					//write the sequence names
					bwTEMP.write(seqName.substring(1)+"\t");
					nGoodTaxa++;
				}
				if(nLociEach[i]>0){
					bwTEMP.write(seq+"\n");
				}
			}
			br.close();
			bw.flush();
			bw.close();

			bwTEMP.write("\n");
			firstGoodLocus=false;

		}
		System.out.println();
		bwX.flush();
		bwX.close();
		bwCharSets.flush();
		bwCharSets.close();
		bwTEMP.flush();
		bwTEMP.close();
		

		BufferedWriter bw2 = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "../RAxML"+folderTail+"/"+project+"_ConcatLoci.phylip" )) ));
		bw2.write(nGoodTaxa+"\t"+lastCharPos+"\n");

		BufferedReader brTEMP = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File( "temp" )) ));
		String tempS=brTEMP.readLine();
		while(tempS!=null){
			bw2.write(tempS+"\n");
			tempS=brTEMP.readLine();
		}
		brTEMP.close();
		new File("temp").delete();
		bw2.flush();
		bw2.close();

		BufferedWriter bw9 = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(   "../Results/"+project+"_AlignmentSummary_SequencePresenceInConcatPhylip"+folderTail+".txt" )) ));
		for(int i=0; i<nTaxa; i++){
			for(int loc=1; loc<=nLoci; loc++){
				if(seqMissing[i][loc-1]){bw9.write("0\t");}else{bw9.write("1\t");}
			}
			bw9.write("\n");
		}
		bw9.flush();
		bw9.close();

      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>> "+ioe.getMessage());}
  }

}
