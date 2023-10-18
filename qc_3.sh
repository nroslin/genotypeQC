#!/bin/bash

module load plink/1.90b6.21
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
scriptdir=/hpf/projects/arnold/users/nroslin/Scripts/Genotypes/qc



###########
# creating a new fam file that includes inferred sex when real sex info is not avail 

awk '{$3==0?sex=$4:sex=$3} NR>1 {print $1,$2,sex}'    ${prefix}_inferredSex.txt > ${tmpfile}_update_sex 
awk '$3==2 {print $1,$2}' ${tmpfile}_update_sex > ${tmpfile}_females 

plink --memory 8000  --bfile $dir/$prefix --update-sex ${tmpfile}_update_sex  --make-bed --out ${tmpfile}_update_sex --set-hh-missing 
# 

# outlier file are non-EUR samples (more precisely, outliers wrt pca)
#cat ${prefix}.remove.txt ${prefix}.1kgpca.outliers.txt| sort -u > ${tmpfile}.remove 
# AS OF v0.1.5 I AM USING CLOSEST 1KG ANCESTRY TO DEFINE GROUPS 


n=`wc -l $prefix.1kgpca.closestAncestry.txt | awk '{print $1-1}'`
groups=`awk 'NR>1 {print $3}' $prefix.1kgpca.closestAncestry.txt | sort -u` 

for g in $groups; do  

 nanc=`awk '$3=="'$g'" {print}' $prefix.1kgpca.closestAncestry.txt | wc -l | awk '{print $1}'` 
 
 if [ $nanc -gt 49 ]; then 
   awk  '$3!="'$g'" {print $1,$2}' $prefix.1kgpca.closestAncestry.txt  > ${tmpfile}.remove 
   plink --memory 8000  --bfile $dir/$prefix  --remove ${tmpfile}.remove --hardy --out ${prefix}_hardy_${g} 

  # removing chrX  
  awk '$1!=23 {print}' ${prefix}_hardy_${g}.hwe > ${tmpfile}.hwe
  \mv ${tmpfile}.hwe ${prefix}_hardy_${g}.hwe 

  # chrX-specific - females only 

  plink --memory 8000  --bfile ${tmpfile}_update_sex  --remove ${tmpfile}.remove --keep ${tmpfile}_females --hardy --out ${tmpfile} --chr 23 

# BUG IN v0.2.4 THE HEADER WAS APPENDED AS WELL SO R SCRIPT WOULD NOT WORK 
  cat ${tmpfile}.hwe | awk 'NR>1 {print}' >> ${prefix}_hardy_${g}.hwe

  R --no-save --args ${prefix}_hardy_${g}.hwe ${prefix}.exclude.txt  < ${scriptdir}/qc_hwe.r 

 fi 
done 


#new to v0.6.0
#add any non-TAG sample to the list of samples to remove (should be controls)
awk '$2!~/^TAG/ {print $1,$2, "CONTROL"}' $dir/$prefix.fam >> $prefix.remove.txt



\rm -f ${tmpfile}* 
\rm -f ${prefix}*.nosex
\rm -f ${prefix}*.nof
\rm -f ${prefix}*.hh






















#############################################################################
#   Copyright 2019 Mathieu Lemire
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#############################################################################
