#!/bin/sh


# Created 31 May 2023.
# Last modified:  21 Jun 2023

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
label=$4   #file with labels for batches (one label per line)
scriptdir=$5

#samples and SNPs to remove
remove=${prefix}.remove.txt
exclude=${prefix}.exclude.txt

#do association test, trend test only, minimum cell count of 5
plink --memory 8000 --allow-no-sex --bfile $dir/$prefix --remove $remove --exclude $exclude --model trend-only --cell 5 --pheno $batch --all-pheno --chr 1-22 --out ${prefix}.batch
  #makes $prefix.batch.PHENNAME.model

#get list of output files
ls *.model > model.tmp$$
nfiles=`wc -l model.tmp$$ | awk '{print $1}'`
i=1
while [ $i -le $nfiles ]
do
  file=`head -$i model.tmp$$ | tail -1`
  name=`head -$i $label | tail -1`

  R --no-save --args $file $name < $scriptdir/qc_batch_qq.r
  i=`expr $i + 1`
done

rm -f ${prefix}.batch.hh ${prefix}.batch.nosex
rm model.tmp$$


#maybe also do PCA, colouring by batch
#what about ancestry and relatedness?
