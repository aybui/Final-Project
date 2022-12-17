# Final-Project
# Phylogenetic analyses of the Chalcid family Encyrtidae using AHE data

#### Folders "AA", "NT12", "NT123", "NTALL" contain the final results and the files for the phylogenetic trees
#### Folder "Ashley" has the code and the reference files

#### "T28.txt" contains the I#####s associated with each specimen used and their species
#### "allChalRefs.txt" are the list of reference transcriptomes that we will use
#### "ReferenceList.txt" contains the I#####s of the refernce chalcids that we'll actually use
#### the I##### folders contain a "I#####_conSeqs.fasta" file (a fasta file which has the consensus sequences) and a "I#####_nMapped.txt" file

#### "T28.AHE(edited).sh" is where the main script is. Most changes and the code I wrote in this script are to adjust for the inclusion of the NT12 (codon w/o wobble base) and NT123 (coding regions only) and changing the script from creating the tree using RAsML to IQTREE
#### "T28.AHE(edited).sh" utilizes scripts from the "Code" 

#### The other script I edited and changes was "Final-Project/Ashley/Code/TranslateHomologs3.java"
#### Was edited and written to include and produce the necessay files to analyze NT123 and NT12 

#### Created new "Align_AllViaMAFFT....sh" scripts for NT12 and NT123 analyses
