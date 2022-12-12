import java.io.*;
import java.util.*;

public class EvalInformativeSitesInPhylipMatrix{
  public static void main(String[] args){
      try{
      
        String infile=args[0];
        
  	BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(infile) ) ));
  	BufferedWriter bwDiff = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(infile+".pairDiff.txt") ) ));
  	BufferedWriter bwSame = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(infile+".pairSame.txt") ) ));

	String tempS=br.readLine();

	int ntaxa = Integer.parseInt(tempS.split("\t")[0]);
	int nchar = Integer.parseInt(tempS.split("\t")[1]);

	char matrix[][] = new char[ntaxa][nchar];

	int pairwiseDiff[][] = new int[ntaxa][ntaxa];
	int pairwiseSame[][] = new int[ntaxa][ntaxa];
	
//	tempS=br.readLine();	//blank space
	tempS=br.readLine();	//first sequence

	int currSite[] = new int[ntaxa];
	boolean firstGene=true;
	int geneCount=0;

	while(tempS!=null && tempS.length()>0){
		geneCount++;
		for(int taxon=0; taxon<ntaxa; taxon++){
			if(firstGene){
				tempS=tempS.split("\t")[1];
			}
			for(int site=0; site<tempS.length(); site++){
				matrix[taxon][currSite[taxon]]=tempS.charAt(site);
				currSite[taxon]++;
			}
			tempS=br.readLine(); // read next sequence (or space if at end of gene)
		}
		firstGene=false;
			tempS=br.readLine(); //first or seq of next gene or last space
	}
	System.out.println("nGenes="+geneCount);
	int nVariable=0;
	int nInformative=0;
	int nA=0;
	int nT=0;
	int nC=0;
	int nG=0;
	int nATCG=0;
	
	int nMissingChars=0;
	int nPartialMissingChars=0;
	int nTotalChars=0;
	
	for(int site=0; site<nchar; site++){
		nA=nT=nC=nG=nATCG=0;
		for(int taxon=0; taxon<ntaxa; taxon++){
//if(taxon==71){continue;}	//use this trick to show statistics without some taxa
			nTotalChars++;
			
			switch(matrix[taxon][site]){
				case 'A': nA++; nATCG++; break;
				case 'T': nT++; nATCG++; break;
				case 'C': nC++; nATCG++; break;
				case 'G': nG++; nATCG++; break;
				case 'N': nMissingChars++; break;
				case '-': nMissingChars++; break;
				default: nPartialMissingChars++;
			}
			for(int taxonB=0; taxonB<ntaxa; taxonB++){
				if((matrix[taxon][site]=='A' || matrix[taxon][site]=='T' || matrix[taxon][site]=='C' || matrix[taxon][site]=='G') && (matrix[taxonB][site]=='A' || matrix[taxonB][site]=='T' || matrix[taxonB][site]=='C' || matrix[taxonB][site]=='G')){
					if(matrix[taxon][site]==matrix[taxonB][site]){
						pairwiseDiff[taxon][taxonB]++;
					}else{
						pairwiseSame[taxon][taxonB]++;					
					}
				}

			}
			
		}
		//check if variable site
		if(nA<nATCG && nT<nATCG && nC<nATCG && nG<nATCG){	//if none of the base counts equal the total then more than one must be greater than zero ** this ignores missing data and N's
			nVariable++;
		}
		if(nA>1 && nA<(nATCG-1) || nT>1 && nT<(nATCG-1) || nC>1 && nC<(nATCG-1) || nG>1 && nG<(nATCG-1)){	//informative if any one of the base counts is between 1 and N-1
			nInformative++;
		}
	}

	System.out.println("Number of taxa = "+ntaxa);
	System.out.println("Number of sites = "+nchar);
	System.out.println("Number of variable sites = "+nVariable);
	System.out.println("Number of informative sites = "+nInformative);
	System.out.println("Number of characters (total) = "+nTotalChars);
	System.out.println("% Missing characters (N's and -'s considered only) ="+(100*nMissingChars/(double)(nTotalChars)));
	System.out.println("% Missing characters (all considered) ="+(100*(nMissingChars+nPartialMissingChars)/(double)(nTotalChars)));
	
	for(int i=0; i<ntaxa; i++){
		for(int j=0; j<ntaxa; j++){
			bwDiff.write(pairwiseDiff[i][j]+"\t");
			bwSame.write(pairwiseSame[i][j]+"\t");
		}
		bwDiff.write("\n");
		bwSame.write("\n");
	}
	bwDiff.flush();
	bwSame.flush();
	bwDiff.close();
	bwSame.close();
	

      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}
  }


}
