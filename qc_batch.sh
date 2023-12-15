#!/bin/sh


# Created 31 May 2023.
# Last modified:  07 Dec 2023

# Optional QC step:  looking for batch/plate effects.  Testing a sliding
# window of plates.  [If examine one vs. the rest, then the test is 
# overpowered; ie, the others combined is a very large sample size.]
#
# Note that if people of different ancestries are in different batches, it 
# may look like a batch effect, so only test in a single ancestry.

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
label=$4   #file with labels for batches (one label per line)
keep=$5   #IDs to keep, from a consistent ancetry
scriptdir=$6

#samples and SNPs to remove
remove=${prefix}.remove.txt
exclude=${prefix}.exclude.txt

#approximate PCs, use the first 4 as covariates
pca=${prefix}.qc_pca.eigenvec

#do some LD pruning
plink  --memory 8000 --bfile $dir/$prefix --remove $remove --exclude $exclude --keep $keep --indep-pairwise 1000 300 0.3 --out tmp$$.ld
plink  --memory 8000 --bfile $dir/$prefix --remove $remove --exclude $exclude --keep $keep --extract tmp$$.ld.prune.in --make-bed --out tmp$$.pruned

#do association test, trend test only, minimum cell count of 5
plink --memory 8000 --allow-no-sex --bfile tmp$$.pruned --logistic --cell 5 --pheno $batch --all-pheno --chr 1-22 --covar $pca --covar-number 1-4 --out ${prefix}.batch
  #makes $prefix.batch.PHENNAME.assoc.logistic
#--covar $pca --covar-number 1-4

#get list of output files
#problem:  if just list files, will be unix sorted, but $label file
#is numerically sorted
#will assume there is no header in the batch file, so that PLINK output files
#will always be P1, P2, etc.
nfiles=`wc -l $label | awk '{print $1}'`
i=1
while [ $i -le $nfiles ]
do
  file=${prefix}.batch.P$i.assoc.logistic
  name=`head -$i $label | tail -1`

  R --no-save --args $file $name < $scriptdir/qc_batch_qq.r
  i=`expr $i + 1`
done

rm -f ${prefix}.batch.hh ${prefix}.batch.nosex
rm -f model.tmp$$
rm -f tmp$$.ld* tmp$$.pruned*

if [ ! -d Batch ]
then
  mkdir Batch
fi
mv *.png Batch
mv *.logistic Batch


#maybe also do PCA, colouring by batch?
#what about and relatedness?
