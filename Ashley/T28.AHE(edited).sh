
#Go to the directory with all the code files:

cd /Users/heraty/Desktop/Ashley/Code
mkdir ../Results

#IDENTIFY CONSENSUS SEQUENCES PASSING COVERAGE FILTER
#Rscript IdentifyGoodSeqsViaReadsMappedD.r $BEGIND $ENDIND 6 100
##Generates minreads.txt file which is the miniumum read file, generates the cutoff for filtering on each exon/locus
Rscript IdentifyGoodSeqsViaReadsMappedD.r 8134 27447 6 100

##Junxia Notes: compile ###.java code to get ###.class before run the java scripts
#javac GatherALLConSeqsWithOKCoverage3.java



#USE HOMOLOG IDS TO OBTAIN GOOD HOMOLOG SEQUENCES TO PUT IN HOMOLOGS FOLDER
#java -Xmx4G GatherALLConSeqsWithOKCoverage3 ../${TAXONSET}.txt $NUMLOCI minReads.txt
java -Xmx4G GatherALLConSeqsWithOKCoverage3 ../T28.txt 1121 minReads.txt
  ## Generates homolog folder w/ consensus seq from all the taxa together for exons that passed our threshold
  ## Removes Contamination identified by Junxia



### ADD REFERENCE TRANSCRIPTOMES TO HOMOLOGS (MANUAL PROCESS)
	# For example, for refs I20942.seqs, I20943.seqs, I20945.seqs, I20946.seqs:
mkdir -p ../HomologsMultipleRefs
for i in {1..1121}; do cp ../Homologs/T28_L$i.fasta ../HomologsMultipleRefs/T28_L$i.fasta && for j in {42,43,49,51,52,54,55,64,76,79,86,87,88,91};
do echo -e ">I209${j}_REFERENCE_Copy1" >> ../HomologsMultipleRefs/T28_L$i.fasta && cat ../References/I209${j}.seqs | head -n $i | tail -n 1
>> ../HomologsMultipleRefs/T28_L$i.fasta; done; done

###Add one reference: I20945.seqs
#mkdir -p ../HomologsOneRef
#for i in {1..1121}; do cp ../Homologs/T10_L$i.fasta ../HomologsOneRef/T10_L$i.fasta && echo -e ">I20945_REFERENCE_Copy1" >> ../HomologsOneRef/T10_L$i.fasta && cat ../References/I20945.seqs | head -n $i | tail -n 1 >> ../HomologsOneRef/T10_L$i.fasta; done
  ##Reference folders has all 49 targeted transciptomes from the lemmons. Each line in the lemmons files corrospond to exon #
  ##adds the reference seq to each of the homolog fasta files and adds sequence headers:


# TRANSLATE TO AA ( CAN USE TRANSLATEHOMOLOGS2 IF ONLY INCLUDING ONE REFERENCE )
	# TranslateHomologs checks for refs in this list using Java's String.startsWith(">"+refID)
java TranslateHomologs3 ../HomologsMultipleRefs/T28_L ../ReferenceList.txt 1121

#Above script generate following folders:
#HomologsAA    (The translated Amino acid sequences)
#HomologsNT12  (Only include first 2 codon posistions; in nucleotides )
#HomologsNT123 (Only the coding region for nucleotide sequences; in nucleotides)
#HomologsNTALL (Both coding and non coding regions; in nucleotides)

#java TranslateHomologs2 ../HomologsOneRef/T10_L I20945 1121
#Works!


### DISTANCES & ORTHOLOGY (Junxia Notes: Do it for all five options!!! This may take a while.)

java GetPairwiseDistanceMeasuresC T15 1121 20 ''
#sometimes get error: Computing distances for locus 459Exception in thread "main" java.lang.OutOfMemoryError: GC overhead limit exceeded
	#at java.util.Arrays.copyOfRange(Arrays.java:3664)
	#at java.lang.String.<init>(String.java:207)
	#at java.lang.StringBuilder.toString(StringBuilder.java:407)
	#at GetPairwiseDistanceMeasuresC.main(GetPairwiseDistanceMeasuresC.java:77)

#Fix: adding in the command line "-XX:-UseGCOverheadLimit" fix the problem:
java -Xmx6144M -XX:-UseGCOverheadLimit GetPairwiseDistanceMeasuresC T28 1121 20 ''
mkdir -p ../OrthoSets
Rscript leanMDS.r T28 1 1121 ''

java GetPairwiseDistanceMeasuresC T28 1121 20 AA
mkdir -p ../OrthoSetsAA
Rscript leanMDS.r T28 1 1121 AA


java GetPairwiseDistanceMeasuresC T28 1121 20 NT12
mkdir -p ../OrthoSetsNT12
Rscript leanMDS.r T28 1 1121 NT12

java GetPairwiseDistanceMeasuresC T28 1121 20 NT123
mkdir -p ../OrthoSetsNT123
Rscript leanMDS.r T28 1 1121 NT123

java GetPairwiseDistanceMeasuresC T28 1121 20 NTALL
mkdir -p ../OrthoSetsNTALL
Rscript leanMDS.r T28 1 1121 NTALL

#ignore the warning message such as:
#1: In finalOrthoIDs == potentialGroup :
#  longer object length is not a multiple of shorter object length


for type in {'',AA,NT12,NT123,NTALL}; do
	java GetPairwiseDistanceMeasuresC T10 1121 20 $type
	mkdir -p ../OrthoSets$type
	Rscript leanMDS.r T10 1 1121 $type > /dev/null &
#done

#SORT THE SEQUENCES INTO ORTHOLOGOUS SETS AND SETUP ALIGNMENTS
	# prelim testing showed that NT12 was best ortho type.  should default to conSeqs (empty string "") as homo type.
	# outputs to Prealignments${HOMOTYPE}${FOLDERTAIL}
	# writes alignment scripts as Align_AllViaMAFFT${HOMOTYPE}${FOLDERTAIL}_conSeqs_*.sh
#java SortSequencesViaOrthoSetsE $TAXONSET ../${TAXONSET}_PloidySpecificationFile_${PLOIDYTAIL}.txt $MINSEQS $NUMTHREADS $FOLDERTAIL $ORTHOTYPE $HOMOTYPE
#Use NT12 distance and pull out NTALL (coding and non-coding regions of seq. that could be translated) for alignment
java SortSequencesViaOrthoSetsE T28 ../T28_PloidySpecificationFile.txt 14 4 _N12N123 NT12 NT123
java SortSequencesViaOrthoSetsE T28 ../T28_PloidySpecificationFile.txt 14 4 _N12NAll NT12 NTAll

#935 loci that have art least 7 sequences


bash Align_AllViaMAFFTNTALL_N12NALL_1.sh
bash Align_AllViaMAFFTNTALL_N12NALL_2.sh
bash Align_AllViaMAFFTNTALL_N12NALL_3.sh
bash Align_AllViaMAFFTNTALL_N12NALL_4.sh

mv ../AlignmentsNTALL_N12NALL ../Alignments_N12NTALL
java TrimAndMaskRawAlignments3 T28 850 28 20 10 0.3 11 ../T28_LociAndTaxaToRemove.txt _N12NTALL
#java SetupRAxML3 TAXSET NEWNLOCI RAXMLJOBS RAXMLJOBS BOOTS BOOTS OUTGROUPTAXON1,OUTGROUPTAXON2 FOLDERTAIL
java SetupRAxML3 T28 850 6 6 100 100 _N12NTALL


##Folder "RAxML_N12NTALL" for gene alignments and concatenated matrix plus partition file for raxml analyses.



#Could repeat but use  NT12 distance and NT123 (only coding region for alignment)
java SortSequencesViaOrthoSetsE T15 ../T15_PloidySpecificationFile.txt 7 4 _N12N123 NT12 NT123

bash Align_AllViaMAFFTNT123_N12N123_1.sh
bash Align_AllViaMAFFTNT123_N12N123_2.sh
bash Align_AllViaMAFFTNT123_N12N123_3.sh
bash Align_AllViaMAFFTNT123_N12N123_4.sh


#ALIGNMENT TRIMMING
#java TrimAndMaskRawAlignments3 TAXSET NEWNLOCI NEWNTAXA WINDOWSIZE MINGOODSITES MINPROPSAME MISSINGALLOWED ../TAXSET_LociAndTaxaToRemove.txt FOLDERTAIL
mv ../AlignmentsNT123_N12NT123 ../Alignments_N12N123
java TrimAndMaskRawAlignments3 T96 1013 28 20 10 0.5 24 ../T96_LociAndTaxaToRemove.txt _N12N123

#java SetupRAxML3 TAXSET NEWNLOCI RAXMLJOBS RAXMLJOBS BOOTS BOOTS OUTGROUPTAXON1,OUTGROUPTAXON2 FOLDERTAIL
java SetupRAxML3 T28 1013 6 6 100 100 _N12N123

##Folder "RAxML_N12NT123" for gene alignments and concatenated matrix plus partition file for raxml analyses.
