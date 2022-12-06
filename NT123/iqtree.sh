#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:15:00
#SBATCH --output=Ashleys_Builogeny
#SBATCH --mail-user=abui052@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="T28-Encyrtidae"
#SBATCH -p intel

module load iqtree/2.2.1
iqtree2 -s T28_ConcatLoci.phylip -m GTR+I+G -p T28_ConcatLoci.charSets -B 1000
