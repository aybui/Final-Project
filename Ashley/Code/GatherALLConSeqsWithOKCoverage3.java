import java.io.*;
import java.util.*;

////////////////////////////////////
//THIS PROGRAM COLLECTS THE SEQUENCES WITH THE GREATEST COVERAGE (NREADS MAPPED) TO USE IN ALIGNMENTS ACROSS SPECIES
//This version removes sequences with N

public class GatherALLConSeqsWithOKCoverage3{
  public static void main(String[] args){
      try{

		String taxonSetFile = args[0];						//e.g. ../T1.txt   //this file contains the taxon ids for the samples to be included in the alignment
		int nLoci = Integer.parseInt(args[1]);				//e.g. 403
		String minReadsFile = args[2];
		
		new File("../Homologs").mkdir();

		BufferedWriter bw[] = new BufferedWriter[nLoci];	
		for(int loc=0; loc<nLoci; loc++){
			bw[loc] = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../Homologs/"+taxonSetFile.substring(taxonSetFile.lastIndexOf("/")+1,taxonSetFile.lastIndexOf("."))+"_L"+(loc+1)+".fasta") ) ));	//in file
		}
			
		//count the number of taxa included
		BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(taxonSetFile) ) ));	//in file							
		String tempS=br.readLine();
		int nInds=0;
		while(tempS!=null){
			nInds++;
			tempS=br.readLine();
		}
		br.close();

		//read in min reads array from file
		BufferedReader brMinReads = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(minReadsFile) ) ));
		int[] minReads = new int[nInds];
		for(int i = 0; i < nInds; i++)
		{
			String string = brMinReads.readLine();
			System.out.println(string);
			minReads[i] = (int)Double.parseDouble(string);
		}

		//SETUP THE ARRAYS TO HOLD THE RESULTS
		String taxonSet[]=new String[nInds];
		String indSet[]=new String[nInds];
		
		//get the taxon identification
		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(taxonSetFile) ) ));	//in file
		for(int i=0; i<nInds; i++){
			taxonSet[i]=br.readLine();			//e.g. I4010_Cyperaceae_Carex_scoparia
			indSet[i]=taxonSet[i].split("_")[0];
		}
		br.close();

		for(int i=0; i<nInds; i++){
			if(!new File("../"+indSet[i]+"/"+indSet[i]+"_conSeqs.fasta").exists()){
				System.out.println("ERROR: Taxon "+indSet[i]+" folder missing conSeqs file."); return;
			}
			//load up the con seq file
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+indSet[i]+"/"+indSet[i]+"_conSeqs.fasta") ) ));	//in file
			tempS=br.readLine();	//e.g. >L1.1
			
			
			if(indSet[i].startsWith("R")){	//if reference copy all sequences (don't apply coverage threshold)
				while(tempS!=null){
					int loc=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));			
					int sup=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1));			

					bw[loc-1].write(">"+taxonSet[i]+"_Copy"+sup+"\n");	//write header  e.g >I8303_STA_0105_Cyperaceae_Carex_pulicaris_Copy1
					tempS=br.readLine();
					bw[loc-1].write(tempS+"\n");	//write sequence
					tempS=br.readLine();			//next header										
				}			
			}else{
				if(!new File("../"+indSet[i]+"/"+indSet[i]+"_nMapped.txt").exists()){
					System.out.println("ERROR: Taxon "+indSet[i]+" folder missing nMapped.txt file."); continue;
				}
				BufferedReader br2 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+indSet[i]+"/"+indSet[i]+"_nMapped.txt") ) ));	//in file
				String tempS2=br2.readLine();
				while(tempS!=null){
					int nMapped = Integer.parseInt(tempS2.split("\t")[2]);
					int loc=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));	
					int sup=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1));
					if(!tempS2.equals(""+loc+"\t"+sup+"\t"+nMapped)){
						System.out.println("ERROR: conSeq header does not match up with nMapped information."); 
						System.out.println(indSet[i]);
						System.out.println(tempS);
						System.out.println(tempS2);
						return;
					}

					tempS=br.readLine();
					if(nMapped>=minReads[i] && tempS.replace("N","").length()>20){ 	//need to write
						bw[loc-1].write(">"+taxonSet[i]+"_Copy"+sup+"\n");	//write header  e.g >I8303_STA_0105_Cyperaceae_Carex_pulicaris_Copy1
						bw[loc-1].write(tempS+"\n");	//write sequence				
					}else{
					}
					tempS2=br2.readLine();
					tempS=br.readLine();
				}
				br2.close();
			}
			br.close();
		}
		for(int loc=0; loc<nLoci; loc++){
			bw[loc].flush();
			bw[loc].close();
		}

      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}
  }

}
