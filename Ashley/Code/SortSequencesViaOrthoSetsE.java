import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

/*
VERSION E is designed to work with the new orthology pipeline and TranslateHomologs2.java
	-pulls sequences from specified files
	-format must be like homologs files, unaligned fasta
*/

//VERSION B allows the user to utilize conSeqs, any of the phased allele sets. It also deals with the changover to the Taxon set number system (does not use the Project number)
////////////////////////////////////
//THIS PROGRAM TAKES THE ORTHO SET INFORMATION AND CREATES NEW conSeq FILES FOR EACH INDIVIDUAL THAT CAN BE USED IN THE GetPairwiseDistanceMeasures.java program
////////////////////////////////////

public class SortSequencesViaOrthoSetsE{
  public static void main(String[] args){
      try{

		String taxonSet = args[0];							//e.g. T1
		String seqFileSpec = args[1];					//e.g. T1_PloidySpecificationFile.txt
		int minSeqsInOrthoSet = Integer.parseInt(args[2]);	//e.g. 25
		int nAlignmentScripts = Integer.parseInt(args[3]);	//e.g. 8
		String folderTail="";
		if(args.length>4){folderTail=args[4];}			//e.g. _2alleles or _conSeqs or _mixturePloidy
		String orthotype = "";
		if(args.length>5){orthotype=args[5];}
		String homotype = "";
		if(args.length>6){homotype=args[6];}

		//COUNT THE NUMBER OF LOCI (MAXIMUM OBSERVED)
		int nLoci=0;
		for(int loc=1; loc<100000; loc++){
			if(new File("../OrthoSets/"+taxonSet+"_L"+loc+".txt").exists()){nLoci=loc;}
		}
		System.out.println(nLoci+" loci found.");
		if(nLoci==0){return;}
		
		//COUNT THE NUBMER OF TAXA IN THE TAXON SET, AND CHECK TO SEE IF ANY ARE MISSING
		BufferedReader brT = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+taxonSet+".txt") ) ));
		String tempT=brT.readLine();	//e.g. I0463_Squamata_Teiidae_Tupinambis_merianae_FG059
		int nTaxa=0;
		while(tempT!=null){
			String tax=tempT.split("_")[0];
			if(!new File("../"+tempT.split("_")[0]).exists()){System.out.println("WARNING: Taxon "+tempT+" not found!");}
			nTaxa++;
			tempT=brT.readLine();
		}
		brT.close();	
		
		//obtain the taxon names from the taxon set file, need this to fill in missing taxa before alignment step		
		String taxa[] = new String[nTaxa];
		brT = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+taxonSet+".txt") ) ));
		for(int i=0; i<nTaxa; i++){
			taxa[i]=brT.readLine();	//e.g. I2982_CER85.2_Gobiidae_Chlamydogobius_eremius
		}
		brT.close();

		//obtain the ploidy specification (list of sequence files to use)
		String seqFiles[] = new String[nTaxa];
		brT = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(seqFileSpec) ) ));
		int ploidy[] = new int[nTaxa];
		for(int i=0; i<nTaxa; i++){
			seqFiles[i]=brT.readLine().split("\t")[0];	//e.g. I2982_allelesFromPost_2.txt

			BufferedReader brP = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+taxa[i].split("_")[0]+"/"+seqFiles[i] )) ));
			//open up each file and determine the ploidy
			String tempP = brP.readLine();	//first header
			if(tempP.split("\\.").length==2){ploidy[i]=1; continue;}
			while(tempP!=null){
				ploidy[i]=Math.max(ploidy[i],Integer.parseInt(tempP.split(".")[2]));
				tempP=brP.readLine();	//sequence
				tempP=brP.readLine();	//next header
			}
		}
		brT.close();		

		//CREATE THE DIRECTORY STRUCTURE FOR BOTH THE PREALIGNMENTS AND THE ALIGNMENTS
		new File("../Prealignments"+homotype+""+folderTail).mkdir();
		new File("../Alignments"+homotype+""+folderTail).mkdir();

		//Setup conversion file
		BufferedWriter bwKey = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Prealignments"+homotype+""+folderTail+"/LocusConversionKEY.sh") ) ));	//out file
		bwKey.write("OldLocNumb\tNewLocNumb\n");
		
		int nLociNew=0;
		for(int loc=1; loc<=nLoci; loc++){
			System.out.println("L"+loc);
			if(!new File("../OrthoSets"+orthotype+"/"+taxonSet+"_L"+loc+".txt").exists()){continue;}
			BufferedReader br2 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../OrthoSets"+orthotype+"/"+taxonSet+"_L"+loc+".txt") ) ));	//file with ortho set IDS.
			//determine how many orthologs sets there could be for this locus
			String tempS2=br2.readLine();
			int currValue=0;
			int max=0;
			while(tempS2!=null){
				currValue=Integer.parseInt(tempS2);
				if(currValue>max){max=currValue;}
				tempS2=br2.readLine();
			}
			br2.close();

			//go through file again and tally up number of seqs in each set
			int tally[] = new int[max];
			br2 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../OrthoSets"+orthotype+"/"+taxonSet+"_L"+loc+".txt") ) ));	//file with ortho set IDS.
			tempS2=br2.readLine();			
			while(tempS2!=null){
				currValue=Integer.parseInt(tempS2);
				tally[currValue-1]++;
				tempS2=br2.readLine();
			}
			br2.close();

			//READ IN HOMOLOGS AND ORTHOSETS FILES
			BufferedReader brH = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../Homologs"+homotype+"/"+taxonSet+"_L"+loc+".fasta") ) ));
			ArrayList<String> homologs = new ArrayList<String>();
			String homolog = null;
			while( (homolog = brH.readLine()) != null)
			{
				homologs.add(homolog);
			}

			BufferedReader brO = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../OrthoSets"+orthotype+"/"+taxonSet+"_L"+loc+".txt") ) ));	
			ArrayList<Integer> orthoIDs = new ArrayList<Integer>();
			String orthoID = null;
			while( (orthoID = brO.readLine()) != null)
			{
				orthoIDs.add(Integer.parseInt(orthoID));
			}

			//TRANSFER SEQUENCES TO PREALIGNMENTS FILES
			for(int i = 0; i < tally.length; i++)
			{
				if(tally[i] >= minSeqsInOrthoSet)
				{
					nLociNew++;

					BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Prealignments"+homotype+""+folderTail+"/"+taxonSet+"_L"+nLociNew+".fasta") ) ));
					BufferedWriter bwKeys = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Prealignments"+homotype+""+folderTail+"/"+taxonSet+"_L"+nLociNew+"_keys.txt") ) ));

					bwKey.write(loc+"\t"+nLociNew+"\n");

					ArrayList<String> sequencesToWrite = new ArrayList<String>();
					ArrayList<String> keysToWrite = new ArrayList<String>();
					for(int j = 0; j < taxa.length; j++)
					{
						for(int a = 1; a <= ploidy[j]; a++)
						{
							sequencesToWrite.add(">"+taxa[j]+"_seq"+a);
							sequencesToWrite.add("N");
							keysToWrite.add(taxa[j]+"\t"+a+"\tNA");
						}
					}

					for(int j = 0; j < orthoIDs.size(); j++)
					{			
						if(orthoIDs.get(j) == i+1)
						{
							String homologHeader = homologs.get(2*j);
							String homologSequence = homologs.get(2*j+1);
							String individual = homologHeader.split("_")[0];

							for(int k = 0; k < sequencesToWrite.size(); k+=2)
							{
								if(sequencesToWrite.get(k).contains(individual))
								{
									sequencesToWrite.set(k+1, homologSequence);
									keysToWrite.set(k/2, keysToWrite.get(k/2).split("_seq")[0]+"_seq"+Integer.parseInt(homologHeader.substring(homologHeader.lastIndexOf("_Copy")+5)));
								}
							}
						}
					}

					for(String line : sequencesToWrite)
					{
						bw.write(line + "\n");
					}
					for(String line : keysToWrite)
					{
						bwKeys.write(line + "\n");
					}

					bw.flush();
					bw.close();
					bwKeys.flush();
					bwKeys.close();
				}
			}
		}
		bwKey.flush();
		bwKey.close();		

		System.out.println(""+nLociNew+" orthologus loci with >= "+minSeqsInOrthoSet+" taxa were found.");

		BufferedWriter bwX[] = new BufferedWriter[nAlignmentScripts];
		for(int i=1; i<=nAlignmentScripts; i++){
			bwX[i-1] = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("Align_AllViaMAFFT"+homotype+""+folderTail+"_"+i+".sh") ) ));	//out file
			bwX[i-1].write("unset MAFFT_BINARIES\n");
		}
		for(int i=0; i<nLociNew; i++){
			if(new File("../Prealignments"+homotype+""+folderTail+"/"+taxonSet+"_L"+(i+1)+".fasta").exists()){
				bwX[i%nAlignmentScripts].write("mafft --genafpair --maxiterate 1000 --quiet ../Prealignments"+homotype+""+folderTail+"/"+taxonSet+"_L"+(i+1)+".fasta > ../Alignments"+homotype+""+folderTail+"/"+taxonSet+"_L"+(i+1)+".fasta\n");
			}
		}
		for(int i=0; i<nAlignmentScripts; i++){		
			bwX[i].flush();
			bwX[i].close();
		}


      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}
  }

}

