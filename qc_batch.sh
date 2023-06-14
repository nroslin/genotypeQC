#!/bin/sh


# Created 31 May 2023.
# Last modified:  31 May 2023

# Optional QC step:  looking for batch/plate effects.  Testing one batch
# vs. the others combined.
# Note that if people of different ancestries are in different batches, it 
# may look like a batch effect.

# Input data:  plink bundle of cleaned genotypes
#	list of IDs and batch number (plain text file)

module load plink/1.90b6.21
module load R/3.6.1

#location of plink files
dir=$1
prefix=$2

#file with batch information
#right now, this needs to be a file that PLINK can interpret as a phenotype file
batch=$3

#samples and SNPs to remove
remove=${prefix}.remove.txt
exclude=${prefix}.exclude.txt

#do association test, trend test only, minimum cell count of 5
plink --memory 8000 --bfile $dir/$prefix --remove $remove --exclude $exclude --model trend-only --cell 5 --pheno $batch --all-pheno --chr 1-23 --out ${prefix}.batch
  #makes $prefix.batch.PHENNAME.model

#making the qq plots in the script might be tricky because of inconsistent
#PHENNAMEs.  Maybe do a grep on *.model files and feed them into the R script
#individually?

#see TCAG desktop /Users/nroslin/Projects/Association/Paterson2017/11Batch

#maybe also do PCA, colouring by batch
#what about ancestry and relatedness?
