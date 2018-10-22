#!/bin/bash

module load plink/1.07
module load R/3.5.1

alias plink='plink --noweb'

# usage, e.g. ./qc.sh ../called SHA15345
# bfiles are expected to be found there:
dir=$1
# prefix of the bfiles
prefix=$2

# tmp file prefix. All _$$_tmp_* will be deleted
tmpfile=_$$_tmp_ 
# where the scripts called here are located 
scriptdir=/hpf/projects/arnold/users/mlemire/scripts/qc/genotypingarrays



imissrate=0.03
lmissrate=0.03

# outlier file are non-EUR samples (more precisely, outliers wrt pca)
cat ${prefix}.remove.txt ${prefix}.1kgpca.outliers.txt| sort -u > ${tmpfile}.remove 


plink  --bfile $dir/$prefix --remove ${tmpfile}.remove --hardy --out ${prefix}_hardy 

R --no-save --args ${prefix}_hardy.hwe ${prefix}.exclude.txt  < ${scriptdir}/qc_hwe.r 


\rm ${tmpfile}* 
\rm ${prefix}*.nosex
\rm ${prefix}*.nof
\rm ${prefix}*.hh













