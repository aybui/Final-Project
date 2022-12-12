BEGIND=
ENDIND=
TAXONSET=
NUMLOCI=
REFERENCELIST=
PLOIDYTAIL=
FOLDERTAIL=
ORTHOTYPE=
HOMOTYPE=
MINSEQS=
NUMTHREADS=


###Junxia Notes: we don't need to do the "SETUP" step.
# SETUP
# Results
#mkdir Code
#cd Code
#cp '~/Dropbox/Anchor_Bioinformatics/Code/AAOrtho_Code.tar.gz' . && tar -zxvf AAOrtho_Code.tar.gz && mv AAOrtho_Code/* . && rmdir AAOrtho_Code && javac *.java


#Go to the directory with all the code files:

cd /Users/jxzhang/Desktop/NewAHE.Pipeline.training.20190206/Code 
mkdir ../Results

#IDENTIFY CONSENSUS SEQUENCES PASSING COVERAGE FILTER
#Rscript IdentifyGoodSeqsViaReadsMappedD.r $BEGIND(first I####) $ENDIND(last I####) 6 100 Gets stats for read mapped to each exon. Generates minimum read file for filtering from I#### Fasta
Rscript IdentifyGoodSeqsViaReadsMappedD.r 5011 26435 6 100
# Output should look this this
	#r 5011 26435 6 100
	#[1] TRUE
	#Generates /Results/GoodConSeqIDs.txt

##Junxia Notes: compile ###.java code to get ###.class before run the java scripts (When you get modifications from Jesse they come in java format. Need to compile into system create .class)
#javac GatherALLConSeqsWithOKCoverage3.java

#USE HOMOLOG IDS TO OBTAIN GOOD HOMOLOG SEQUENCES TO PUT IN HOMOLOGS FOLDER (Taxon based filtering cutoff. Undercut off, consensus sequence are filtered out based on coverage.)
#java -Xmx4G GatherALLConSeqsWithOKCoverage3 ../${TAXONSET}.txt(Taxon Folder, with I####s and names) $NUMLOCI(Number exon recovered in assembly) minReads.txt (Tells script where to find the minimum reads)
java -Xmx4G GatherALLConSeqsWithOKCoverage3 ../T15.txt 1121 minReads.txt
#Output should look like
	#500
	#500
	#500
	#500
	#500
	#8
	#5
	#33
	#2
	#4
	#151
	#35
	#6
	#8
	#3
#Generates Homologs folder with consensus sequences from all taxa together for exons that passed our taxon specific threshold.



### ADD REFERENCE TRANSCRIPTOMES TO HOMOLOGS (MANUAL PROCESS) (References folder has all 49 targeted transcriptomes from the lemmons. Each line in file corrosponds to exon one, line to to exon 2 etc. Adds reference sequences into each of the homolog files. Make sure to use multiple reference sequences. Select as many reference sequences you can within your group, or closely allied groups.)
	# For example, for refs I20942.seqs, I20943.seqs, I20945.seqs, I20946.seqs:
mkdir -p ../HomologsMultipleRefs
#fori in {}; do cp ../{HOMOLOGS}.fasta {TOHOMOLOGSMULTIPLEREFS}.fasta && for j in {reference_trascriptomes_last_2 digits}; do echo -e ">I209${j}_REFERENCE_Copy1" >> {}
for i in {1..1121}; do cp ../Homologs/T15_L$i.fasta ../HomologsMultipleRefs/T15_L$i.fasta && for j in {42,43,45,46}; do echo -e ">I209${j}_REFERENCE_Copy1" >> ../HomologsMultipleRefs/T15_L$i.fasta && cat ../References/I209${j}.seqs | head -n $i | tail -n 1 >> ../HomologsMultipleRefs/T15_L$i.fasta; done; done

###Add one reference: I20945.seqs
#mkdir -p ../HomologsOneRef
#for i in {1..1121}; do cp ../Homologs/T10_L$i.fasta ../HomologsOneRef/T10_L$i.fasta && echo -e ">I20945_REFERENCE_Copy1" >> ../HomologsOneRef/T10_L$i.fasta && cat ../References/I20945.seqs | head -n $i | tail -n 1 >> ../HomologsOneRef/T10_L$i.fasta; done



# TRANSLATE TO AA ( CAN USE TRANSLATEHOMOLOGS2 IF ONLY INCLUDING ONE REFERENCE ) (ATTEMPT TO IDENTIFY CODING REGIONS AND TRANSLATE TO AA)
	# TranslateHomologs checks for refs in this list using Java's String.startsWith(">"+refID)
java TranslateHomologs3 ../HomologsMultipleRefs/T15_L ../ReferenceList 1121 
#Reference list has text file with the I###s used for this analysis (needs to be updated to reflect the reference sequence that you are using; NUMLOCI(number of exons that are in analysis))
#(TranslateHomolog3 can handle multiple reference sequences. TranslateHomolog2 can only take one reference folder)

#Above script generate following folders:
#HomologsAA    (The translated Amino acid sequences)
#HomologsNT12  (Only include first 2 codon posistions; in nucleotides )
#HomologsNT123 (Only the coding region for nucleotide sequences; in nucleotides)
#HomologsNTALL (Both coding and non coding regions; in nucleotides)

#java TranslateHomologs2 ../HomologsOneRef/T10_L I20945 1121
#Works!


### DISTANCES & ORTHOLOGY (Junxia Notes: Do it for all five options!!! This may take a while.) (Calculate distance using consensus sequence, AA, NT12, NT123, and NTALL)

java GetPairwiseDistanceMeasuresC T15 1121 20 '' #(Taxonset, Number of exons, sliding window of 20bp or AA, ''(taking consensus sequences))
#sometimes get error: Computing distances for locus 459Exception in thread "main" java.lang.OutOfMemoryError: GC overhead limit exceeded
	#at java.util.Arrays.copyOfRange(Arrays.java:3664)
	#at java.lang.String.<init>(String.java:207)
	#at java.lang.StringBuilder.toString(StringBuilder.java:407)
	#at GetPairwiseDistanceMeasuresC.main(GetPairwiseDistanceMeasuresC.java:77)

#Fix: adding in the command line "-XX:-UseGCOverheadLimit" fix the problem:
java -Xmx6144M -XX:-UseGCOverheadLimit GetPairwiseDistanceMeasuresC T15 1121 20 '' (#Only run when you get the above error)

mkdir -p ../OrthoSets
Rscript leanMDS.r T15 1 1121 ''

#(Run every dataset for the exons)

java GetPairwiseDistanceMeasuresC T15 1121 20 AA
mkdir -p ../OrthoSetsAA
Rscript leanMDS.r T15 1 1121 AA

java GetPairwiseDistanceMeasuresC T15 1121 20 NT12
mkdir -p ../OrthoSetsNT12
Rscript leanMDS.r T15 1 1121 NT12

java GetPairwiseDistanceMeasuresC T15 1121 20 NT123
mkdir -p ../OrthoSetsNT123
Rscript leanMDS.r T15 1 1121 NT123

java GetPairwiseDistanceMeasuresC T15 1121 20 NTALL
mkdir -p ../OrthoSetsNTALL
Rscript leanMDS.r T15 1 1121 NTALL

#ignore the warning message such as:
#1: In finalOrthoIDs == potentialGroup :
#  longer object length is not a multiple of shorter object length


#for type in {'',AA,NT12,NT123,NTALL}; do
	#java GetPairwiseDistanceMeasuresC T10 1121 20 $type
	#mkdir -p ../OrthoSets$type
	#Rscript leanMDS.r T10 1 1121 $type > /dev/null &
#done

#SORT THE SEQUENCES INTO ORTHOLOGOUS SETS AND SETUP ALIGNMENTS
# prelim testing showed that NT12 was best ortho type.  should default to conSeqs (empty string "") as homo type.
	# outputs to Prealignments${HOMOTYPE}${FOLDERTAIL}
	# writes alignment scripts as Align_AllViaMAFFT${HOMOTYPE}${FOLDERTAIL}_conSeqs_*.sh
#java SortSequencesViaOrthoSetsE $TAXONSET ../${TAXONSET}_PloidySpecificationFile_${PLOIDYTAIL}.txt $MINSEQS $NUMTHREADS $FOLDERTAIL $ORTHOTYPE $HOMOTYPE
#Could repeat but use  NT12 distance and NT123 (only coding region for alignment) 
java SortSequencesViaOrthoSetsE T15 ../T15_PloidySpecificationFile.txt 7 4 _N12NT123 NT12 NT123

bash Align_AllViaMAFFTNT123_N12NT123_1.sh
bash Align_AllViaMAFFTNT123_N12NT123_2.sh
bash Align_AllViaMAFFTNT123_N12NT123_3.sh
bash Align_AllViaMAFFTNT123_N12NT123_4.sh


#ALIGNMENT TRIMMING
#java TrimAndMaskRawAlignments3 TAXSET(Taxon Set.txt) NEWNLOCI(exon number in alignment) NEWNTAXA(Number of taxa after filtering) WINDOWSIZE MINGOODSITES (for a 20 base pair window, if X are determined as good, they are kept) MINPROPSAME (if X percent of sites are good they are kept) MISSINGALLOWED (Number of taxa allowed to be missing to trim certain sites) ../TAXSET_LociAndTaxaToRemove.txt FOLDERTAIL
mv ../AlignmentsNT123_N12NT123 ../Alignments_N12NT123
java TrimAndMaskRawAlignments3 T15 935 15 20 10 0.5 9 ../T15_LociAndTaxaToRemove.txt _N12NT123 #Play with the parameters to get the most good data possible.

#java SetupRAxML3 TAXSET NEWNLOCI RAXMLJOBS RAXMLJOBS BOOTS BOOTS OUTGROUPTAXON1,OUTGROUPTAXON2 FOLDERTAIL
java SetupRAxML3 T15 935 6 6 100 100 _N12NT123

##Folder "RAxML_N12NT123" for gene alignments and concatenated matrix plus partition file for raxml analyses.

	
#Use NT12 distance and pull out NTALL (coding and non-coding regions of seq. that could be translated) for alignment
java SortSequencesViaOrthoSetsE T15 ../T15_PloidySpecificationFile.txt 7 4 _N12NTALL NT12 NTALL #(Sorts the orthologs into each sequence file for alignment. Feeds in ploidy, MINSEQS (keep only orthologs with at least this number of taxa) 4(How many parrallel alignments do you want to do. Smaller data set you might want to use NT123 for bigger NT12 to avoid noise NTALL for both coding and non-coding regions NT123 for only coding regions))

#935 loci that have are least 7 sequences

#MAFFT ALIGNMENTS (this will take a while)
bash Align_AllViaMAFFTNTALL_N12NTALL_1.sh
bash Align_AllViaMAFFTNTALL_N12NTALL_2.sh
bash Align_AllViaMAFFTNTALL_N12NTALL_3.sh
bash Align_AllViaMAFFTNTALL_N12NTALL_4.sh


#ALIGNMENT TRIMMING
#java TrimAndMaskRawAlignments3 TAXSET NEWNLOCI NEWNTAXA WINDOWSIZE MINGOODSITES MINPROPSAME MISSINGALLOWED ../TAXSET_LociAndTaxaToRemove.txt FOLDERTAIL

#Before run alignment trimming need to change alignment folder tail: 
mv ../AlignmentsNTALL_N12NTALL ../Alignments_N12NTALL
java TrimAndMaskRawAlignments3 T15 935 15 20 10 0.5 9 ../T15_LociAndTaxaToRemove.txt _N12NTALL

#java SetupRAxML3 TAXSET NEWNLOCI RAXMLJOBS RAXMLJOBS BOOTS BOOTS OUTGROUPTAXON1,OUTGROUPTAXON2 FOLDERTAIL
java SetupRAxML3 T15 935 6 6 100 100 _N12NTALL


##Folder "RAxML_N12NTALL" for gene alignments and concatenated matrix plus partition file for raxml analyses.


