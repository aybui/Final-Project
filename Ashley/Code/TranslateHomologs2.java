
import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
 
//java TranslateHomologs2 ../HomologsToTestOrthology/T000_L I20946 441
public class TranslateHomologs2{
  static int K;
  public static void main(String[] args){
   try{

    //GET INPUT PARAMETERS
    String inFileStem = args[0];            //e.g. Homologs/T461_L
    String refID = args[1];				          //e.g. I20946
    int nLoci = Integer.parseInt(args[2]);  //e.g. 441

    K=15;
	
    for(int loc=1; loc<=nLoci; loc++){

      String infile=inFileStem+loc+".fasta";
      if(!new File(infile).exists()){ continue; }

      //get reference sequence
      BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(infile) ) ));  //in file
      String tempS=br.readLine();
      
      while(tempS!=null && !tempS.startsWith(">"+refID+"_")){
        tempS=br.readLine();
      }
      if(tempS==null){ continue; }  //could not find reference sequence, need to skip this locus...
      String refNT=br.readLine().toUpperCase();
      if(refNT.length() < 3){ continue; }

      // int minStops=999;
      // String minStopAA="";

      // //try each reading frame, count stop codons
      // String refAA1=translate(refNT);
      // int nStops1=0; 
      // for(int x=0; x<refAA1.length(); x++){
      //   if(refAA1.charAt(x)=='*'){nStops1++;}
      // }
      // minStopAA=refAA1;
      // minStops=nStops1;

      // String refAA2=translate(refNT.substring(1));
      // int nStops2=0; 
      // for(int x=0; x<refAA2.length(); x++){
      //   if(refAA2.charAt(x)=='*'){nStops2++;}
      // }
      // if(nStops2<minStops){minStopAA=refAA2; minStops=nStops2;}

      // String refAA3=translate(refNT.substring(2));
      // int nStops3=0; 
      // for(int x=0; x<refAA3.length(); x++){
      //   if(refAA3.charAt(x)=='*'){nStops3++;}
      // }
      // if(nStops3<minStops){minStopAA=refAA3; minStops=nStops3;}

      // System.out.println("L"+loc+"\t"+nStops1+"\t"+nStops2+"\t"+nStops3);
      
      // if(minStops>1){
      //   System.out.println("Warning: Skipping locus "+loc+" because reading frame could not be established for reference sequence (more than one stop codon was found in each reading frame. ");
      //   continue;
      // }        
      // String refAA=minStopAA;

      int numStops = 0;
      String refAA = null;
      for(int i = 0; i < 6; i++)
      {
        if(i<3)
        {
          refAA = translate(refNT.substring(i));
        }
        else
        {
          refAA = translate(getRevComp(refNT).substring(i-3));
        }
        
        numStops = 0;
        for(char c : refAA.toCharArray())
        {
          if(c == '*')
          {
            numStops++;
          }
        }

        if(numStops <= 1)
        {
          break;
        }
      }
      if(numStops > 1)
      {
        System.out.println("Warning: Skipping locus "+loc+" because reading frame could not be established for reference sequence (more than one stop codon was found in each reading frame. ");
      }




      br.close();

      //HASH KMERS FOUND IN REFERENCE SEQUENCE FOR THIS LOCUS
    	HashMap map = new HashMap();
        String kmer="";
        for(int i=0; i<refAA.length()-K+1; i++){
        	kmer=refAA.substring(i,i+K);

        	Integer value = (Integer)map.get(kmer);
        	if(value==null){
        		map.put(kmer, (Integer)i);				//value is position in reference amino acid sequence at which the kmer starts
        	}else{
        		map.put(kmer, (Integer)(-9999));
        	}
        }

        //CHECK THAT THE HOMOLOG FILE EXISTS

    	//DETERMINE THE NUMBER OF SEQUENCES IN THE HOMOLOG FILE
		    //reset the infile
        br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(infile) ) ));  //in file
      	tempS=br.readLine(); //first header
      	int nSeqs=0;
      	while(tempS!=null){
      		nSeqs++;
      		tempS=br.readLine(); //seq
      		tempS=br.readLine(); //next header
      	}
      	br.close();

      	String heads[] = new String[nSeqs];
      	String seqs[] = new String[nSeqs];
      	int lens[] = new int[nSeqs];

    	//LOAD IN ALL OF THE SEQUENCES
        br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(infile) ) ));
        for(int i=0; i<nSeqs; i++){
        	heads[i]=br.readLine();
        	seqs[i]=br.readLine().toUpperCase();
        	lens[i]=seqs[i].length();
        }
        br.close();

        //SETUP THE OUTFILES
        new File("../HomologsAA/").mkdir();
        new File("../HomologsNT12/").mkdir();
        new File("../HomologsNT123/").mkdir();
        new File("../HomologsNTALL/").mkdir();

        String outFileStem = inFileStem.split("/")[inFileStem.split("/").length-1];
        BufferedWriter bwAA = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../HomologsAA/"+outFileStem+loc+".fasta") ) ));  //
        BufferedWriter bwNT12 = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../HomologsNT12/"+outFileStem+loc+".fasta") ) ));  //
        BufferedWriter bwNT123 = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../HomologsNT123/"+outFileStem+loc+".fasta") ) ));  //
        BufferedWriter bwNTALL = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File("../HomologsNTALL/"+outFileStem+loc+".fasta") ) ));  //

        //DETERMINE READING FRAME OF EACH POSITION IN EACH HOMOLOG
        for(int i=0; i<nSeqs; i++){

          if(seqs[i].length() < 3)
          {
            continue;
          }

        	int hits[] = new int[lens[i]];
        	int frame[] = new int[lens[i]];
        	int refPos[] = new int[lens[i]];

        	String tempAA="";

          String homologAA = "";
          String homologNT12 = "";
          String homologNT123 = "";
          String homologNTALL = "";

          int currNAAs=0;
          int bestNAAs=0;

        	//try all three reading frames (assume forward orientation)
        	for(int f=0; f<6; f++){
        		if(f<3){	//forward reading frames
	        		tempAA=translate(seqs[i].substring(f));
	        	}else{		//reverse reading frames
	        		tempAA=translate(getRevComp(seqs[i]).substring(f-3));
	        	}

	        	int lenAA=tempAA.length();

	        	//boolean inFrame[]= new boolean[lenAA];	//default false

	        	Arrays.fill(refPos,-1);
	        	//position in translated amino acid sequence
        		for(int j=0; j<lenAA-K+1; j++){
        			kmer=tempAA.substring(j,j+K);
        			Integer value = (Integer)map.get(kmer);
        			if(value!=null){
        				for(int x=0; x<K; x++){
        					//inFrame[j+x]=true;
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

            //count the number of valid amino acids that were translated for this reading frame...
            currNAAs=0;
            for(int x=0; x<lenAA; x++){
              if(refPos[x]>-1){currNAAs++;}
            }

            if(currNAAs>bestNAAs){  //found better reading frame...store sequences
              bestNAAs=currNAAs;
              //isolate the valid amino acid sequence and corresponding nucleotide sequences
              homologAA = "";
              homologNT123 = "";
              homologNT12 = "";
              homologNTALL = seqs[i];
              for(int x=0; x<lenAA; x++){
                if(refPos[x]>-1){
                  homologAA+=tempAA.charAt(x);
                  homologNT12+=seqs[i].substring((f%3)+x*3,(f%3)+x*3+2);                  
                  homologNT123+=seqs[i].substring((f%3)+x*3,(f%3)+x*3+3);
                }
              }
            }
        	}
        	if(bestNAAs>0){  //we were successful in translating the sequence into one of the reading frames...

            bwAA.write(heads[i]+"\n"+homologAA+"\n");
            bwNT12.write(heads[i]+"\n"+homologNT12+"\n");
            bwNT123.write(heads[i]+"\n"+homologNT123+"\n");
            bwNTALL.write(heads[i]+"\n"+homologNTALL+"\n");
          }
        }

        bwAA.flush(); bwNT123.flush(); bwNT12.flush(); bwNTALL.flush();
        bwAA.close(); bwNT123.close(); bwNT12.close(); bwNTALL.close();

    }

   }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}

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

  static String getRevComp(String seq){

  	String revComp="";
  	int len=seq.length();
  	for(int i=0; i<len; i++){
  		revComp+=seq.charAt(len-(i+1));
  	}
  	if(seq.length()!=revComp.length()){System.out.println("ERROR 1 IN getRevCOMP!!!");}
  
  	revComp=revComp.replace('A','Z');
  	revComp=revComp.replace('T','A');
  	revComp=revComp.replace('Z','T'); 
  	
  	revComp=revComp.replace('C','Q');
  	revComp=revComp.replace('G','C');
  	revComp=revComp.replace('Q','G');
 
  	//REPLACE ALL OTHER BASES WITH 'N'
  	for(int i=0; i<len; i++){
  		char c=revComp.charAt(i);
  		if(c!='-' && c!='A' && c!='T' && c!='C' && c!='G'){
  			if(len>i+1){	revComp=revComp.substring(0,i)+'N'+revComp.substring(i+1);
  			}else{		revComp=revComp.substring(0,i)+'N'; }
		  }
  	}  	
  	if(revComp.length()!=len){System.out.println("ERROR 2 IN getRevCOMP!!!");}
  		
	return revComp;

  }

}


