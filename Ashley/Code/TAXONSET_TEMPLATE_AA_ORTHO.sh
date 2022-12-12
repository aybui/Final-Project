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


# SETUP
mkdir Results
mkdir Code
cd Code
cp '~/Dropbox/Anchor_Bioinformatics/Code/AAOrtho_Code.tar.gz' . && tar -zxvf AAOrtho_Code.tar.gz && mv AAOrtho_Code/* . && rmdir AAOrtho_Code && javac *.java

#IDENTIFY CONSENSUS SEQUENCES PASSING COVERAGE FILTER
Rscript IdentifyGoodSeqsViaReadsMappedD.r $BEGIND $ENDIND 6 100

#USE HOMOLOG IDS TO OBTAIN GOOD HOMOLOG SEQUENCES TO PUT IN HOMOLOGS FOLDER
java -Xmx4G GatherALLConSeqsWithOKCoverage3 ../${TAXONSET}.txt $NUMLOCI minReads.txt

### ADD REFERENCE TRANSCRIPTOMES TO HOMOLOGS (MANUAL PROCESS)
	# For example, for Bond millipede refs myrTra1 - myrTra22:
		#mkdir -p ../HomologsMultipleRefs
		#for i in {1..22}; do cp ../Homologs/T490_L$i.fasta ../HomologsMultipleRefs/T490_L$i.fasta && for j in {1..22}; do echo -e ">myrTra${j}_REFERENCE_Copy1" >> ../HomologsMultipleRefs/T490_L$i.fasta && cat ../References/myrTra${j}_.seqs | head -n $i | tail -n 1 >> ../HomologsMultipleRefs/T490_L$i.fasta; done; done

# TRANSLATE TO AA ( CAN USE TRANSLATEHOMOLOGS2 IF ONLY INCLUDING ONE REFERENCE )
	# TranslateHomologs checks for refs in this list using Java's String.startsWith(">"+refID)
java TranslateHomologs3 ../HomologsMultipleRefs/T490_L $REFERENCELIST 187

### DISTANCES & ORTHOLOGY
for type in {'',AA,NT12,NT123,NTALL}; do
	java GetPairwiseDistanceMeasuresC $TAXONSET $NUMLOCI 20 $type
	mkdir -p ../OrthoSets$type
	Rscript leanMDS.r $TAXONSET 1 $NUMLOCI $type > /dev/null &
done

#SORT THE SEQUENCES INTO ORTHOLOGOUS SETS AND SETUP ALIGNMENTS
	# prelim testing showed that NT12 was best ortho type.  should default to conSeqs (empty string "") as homo type.
	# outputs to Prealignments${HOMOTYPE}${FOLDERTAIL}
	# writes alignment scripts as Align_AllViaMAFFT${HOMOTYPE}${FOLDERTAIL}_conSeqs_*.sh
java SortSequencesViaOrthoSetsE $TAXONSET ../${TAXONSET}_PloidySpecificationFile_${PLOIDYTAIL}.txt $MINSEQS $NUMTHREADS $FOLDERTAIL $ORTHOTYPE $HOMOTYPE

#ALIGNMENTS
chmod +x Align_AllViaMAFFT*
for i in $(seq 1 $NUMTHREADS); do
	./Align_AllViaMAFFT${HOMOTYPE}${FOLDERTAIL}_$i.sh &
done







### END OF AA ORTHO MODIFICATIONS.  Pasting old script below, will clean it up at some point in the future.








########################## TRIMMING ##########################
#note the following block will be repeated until the following parameters are optimized: MINPROPSAME, MINGOODSITES, MISSINGALLOWED
#java TrimAndMaskRawAlignments3 TAXSET NEWNLOCI NEWNTAXA WINDOWSIZE MINGOODSITES MINPROPSAME MISSINGALLOWED ../TAXSET_LociAndTaxaToRemove.txt FOLDERTAIL
java TrimAndMaskRawAlignments3 T490 NEW.N.LOCI NEW.N.TAXA 20 14 0.5 MISSING.ALLOWED ../T490_LociAndTaxaToRemove.txt _conSeqs

#Rscript PlotAlignmentSummary2.r TAXSET NEWNTAXA NEWNLOCI MISSINGALLOWED FOLDERTAIL
Rscript PlotAlignmentSummary2.r "T490" NEW.N.TAXA NEW.N.LOCI MISSING.ALLOWED "_conSeqs"

#Look at plot Results/TAXSET_AlignmentSummary_nLoci.pdf adjust MISSINGALLOWED so it corrsponds to plateau
#Look at plot Results/TAXSET_AlignmentSummary_MaxTally.pdf adjust MINPROPSAME so it corrsponds to valley
#Look at plot Results/TAXSET_AlignmentSummary_WindowTally.pdf adjust MINGOODSITES so it corrsponds to valley

#java SetupRAxML3 TAXSET NEWNLOCI RAXMLJOBS RAXMLJOBS BOOTS BOOTS OUTGROUPTAXON1,OUTGROUPTAXON2 FOLDERTAIL
java SetupRAxML3 T490 NEW.N.LOCI 6 6 100 100 _conSeqs

#Load RAxML/TAXSET_ConcatLoci.phylip in Geneious and inspect the alignment
#Increasing MINGOODSITES will increase masking.
#Increasing MINPROPSAME will increase trimming.

#Remove any bad loci or individuals in TAXSET_LociAndTaxaToRemove.txt
#Check RAxML_FOLDERTAIL/TAXSET_ConcatLoci.charSets to see which locus number needs to be removed
#Format:
LociNumber	Individual
348			ALL				#Removes locus 348 for all individuals.
ALL			21				#Removes all loci for individual 21. This means the 21st individual in Geneious.

#RUN ALL TRIMMING STEPS AGAIN IF BAD SEQUENCES IDENTIFIED



########################## PHYLOGENETICS ##########################
#SUMMARIZE FINAL ALIGNMENT
java EvalInformativeSitesInPhylipMatrix ../RAxML_conSeqs/T490_ConcatLoci.phylip > ../RAxML_conSeqs/T490_Statistics.txt

#OPEN TO SEE THE STATS
cat ../RAxML_conSeqs/T490_Statistics.txt

#CHANGE DIRECTORIES AND PERMISSIONS
cd ../RAxML_conSeqs
chmod +x *.sh



############## CONCATENATED TREES ##############
./T490_ConcatLoci.sh



################ SEPARATE TREES ################
./T490_SeparateLoci.sh



#OPEN RAxML_bipartitions.TAXSET_ConcatLoci IN FIGTREE AND SAVE AS TAXSET_Tree.pdf



################# ASTRAL TREES #################
#MAKE ASTRAL FOLDER
mkdir ../Astral
cd ../Code
rsync -av --progress /home/alemmon/Dropbox/Anchor_Bioinformatics/Code/lib /home/alemmon/Dropbox/Anchor_Bioinformatics/Code/astral.4.11.2.jar .

#CONCATENATE EACH BESTTREE LOCI INTO ONE FILE
cat ../RAxML_conSeqs/RAxML_bestTree.T490_L* > ../Astral/RAxML_bestTree.T490_ALL.txt
#CONCATENATE EACH BOOTSTRAP LOCI INTO ONE FILE 
ls -d -1 ../RAxML_conSeqs/RAxML_bootstrap.T490_L* > ../Astral/bs-files

#RUN ASTRAL
java -jar astral.4.11.2.jar -i ../Astral/RAxML_bestTree.T490_ALL.txt -b ../Astral/bs-files -t 1 -o ../Astral/T490_Astral.tre 2> ../Astral/log

#OPEN IN FIGTREE AND MAKE PDF TREE
tail -n 1 ../Astral/T490_Astral.tre > ../Astral/AstralTree.tre
