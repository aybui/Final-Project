
cd /Desktop/Final-Project/Ashley/Code
mkdir ../Results

#THIS IDENTIFY CONSENSUS SEQUENCES PASSING COVERAGE FILTER
Rscript IdentifyGoodSeqsViaReadsMappedD.r 8134 27447 6 100

#THIS USES HOMOLOG IDS TO OBTAIN GOOD HOMOLOG SEQUENCES TO PUT IN HOMOLOGS FOLDER
java -Xmx4G GatherALLConSeqsWithOKCoverage3 ../T28.txt 1121 minReads.txt
  ## Generates homolog folder w/ consensus seq from all the taxa together for exons that passed our threshold
  ## Removes Contamination

#THIS ADD REFERENCE TRANSCRIPTOMES TO HOMOLOGS
mkdir -p ../HomologsMultipleRefs
for i in {1..1121}; do cp ../Homologs/T28_L$i.fasta ../HomologsMultipleRefs/T28_L$i.fasta && for j in {42,43,49,51,52,54,55,64,76,79,86,87,88,91};
do echo -e ">I209${j}_REFERENCE_Copy1" >> ../HomologsMultipleRefs/T28_L$i.fasta && cat ../References/I209${j}.seqs | head -n $i | tail -n 1
>> ../HomologsMultipleRefs/T28_L$i.fasta; done; done


#THIS TRANSLATE TO AA
## Edited to include NT123, NT12 and NTALL
java TranslateHomologs3 ../HomologsMultipleRefs/T28_L ../ReferenceList.txt 1121


#THIS GETS THE PAIRWISE DISTANCES & ORTHOLOGY
java GetPairwiseDistanceMeasuresC T15 1121 20 ''

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

for type in {'',AA,NT12,NT123,NTALL}; do
	java GetPairwiseDistanceMeasuresC T10 1121 20 $type
	mkdir -p ../OrthoSets$type
	Rscript leanMDS.r T10 1 1121 $type > /dev/null &

###NTALL
java SortSequencesViaOrthoSetsE T28 ../T28_PloidySpecificationFile.txt 14 4 _N12NAll NT12 NTAll

bash Align_AllViaMAFFTNTALL_N12NALL_1.sh
bash Align_AllViaMAFFTNTALL_N12NALL_2.sh
bash Align_AllViaMAFFTNTALL_N12NALL_3.sh
bash Align_AllViaMAFFTNTALL_N12NALL_4.sh

mv ../AlignmentsNTALL_N12NALL ../Alignments_N12NTALL
java TrimAndMaskRawAlignments3 T28 850 28 20 10 0.3 11 ../T28_LociAndTaxaToRemove.txt _N12NTALL
java SetupRAxML3 T28 850 6 6 100 100 _N12NALL
module load iqtree/2.2.1
iqtree2 -s T28_ConcatLoci.phylip -m GTR+I+G -p T28_ConcatLoci.charSets -B 1000

###NT123
java SortSequencesViaOrthoSetsE T28 ../T28_PloidySpecificationFile.txt 14 4 _N12N123 NT12 NT123

bash Align_AllViaMAFFTNT123_N12N123_1.sh
bash Align_AllViaMAFFTNT123_N12N123_2.sh
bash Align_AllViaMAFFTNT123_N12N123_3.sh
bash Align_AllViaMAFFTNT123_N12N123_4.sh

mv ../AlignmentsNT123_N12NT123 ../Alignments_N12N123
java TrimAndMaskRawAlignments3 T96 1013 28 20 10 0.5 24 ../T96_LociAndTaxaToRemove.txt _N12N123
module load iqtree/2.2.1
java SetupRAxML3 T28 850 6 6 100 100 _N12N123
iqtree2 -s T28_ConcatLoci.phylip -m GTR+I+G -p T28_ConcatLoci.charSets -B 1000

###NT12
java SortSequencesViaOrthoSetsE T28 ../T28_PloidySpecificationFile.txt 14 4 _N12N12 NT12 NT12

bash Align_AllViaMAFFTNT12_N12N12_1.sh
bash Align_AllViaMAFFTNT12_N12N12_2.sh
bash Align_AllViaMAFFTNT12_N12N12_3.sh
bash Align_AllViaMAFFTNT12_N12N12_4.sh

mv ../AlignmentsNTALL_N12NALL ../Alignments_N12N12
java TrimAndMaskRawAlignments3 T28 850 28 20 10 0.3 11 ../T28_LociAndTaxaToRemove.txt _N12N12
java SetupRAxML3 T28 850 6 6 100 100 _N12N12
module load iqtree/2.2.1
iqtree2 -s T28_ConcatLoci.phylip -m GTR+I+G -p T28_ConcatLoci.charSets -B 1000
