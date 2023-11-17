#!/bin/bash

# Created 15 November 2023.
# Last modified:  16 Nov 2023

# Process control samples.  Assume if ID does not start with TAG then is a
# control.  Also assume all controls are the same ID (NA12878) and so should
# have identical genotypes.

module load plink/1.90b6.21


echo 
echo "#####################################"
echo "Calling script $0"
echo "Processing of control samples"
echo

# usage, e.g. ./qc.sh ../called SHA15345
# bfiles are expected to be found there:
dir=$1
# prefix of the bfiles
prefix=$2

# tmp file prefix. All _$$_tmp_* will be deleted when this script completes
tmpfile=_$$_tmp_ 

# where the scripts called here are located 
scriptdir=/hpf/projects/arnold/users/nroslin/Scripts/Genotypes/qc

#summary information will be written into a report
reportfile=${prefix}_QCreport.txt
if [ -e $reportfile ]
then
  echo >> $reportfile
else
  echo > $reportfile
fi
echo "Controls" >> $reportfile
echo "--------" >> $reportfile
echo "Processing of control samples" >> $reportfile

awk '$2!~/^TAG/ {print $1,$2}' $dir/$prefix.fam > $tmpfile.control
ncontrols=`wc -l $tmpfile.control | awk '{print $1}'`

if [ $ncontrols -gt 0 ]
then
  #there are controls
  #try to merge them so that we can get concordance rate
  #extract controls
  plink --bfile $dir/$prefix --exclude $prefix.exclude.txt --keep $tmpfile.control --make-bed --out ${tmpfile}_control

  #rename so they all have the same ID
  mv ${tmpfile}_control.fam ${tmpfile}_old.fam
  awk '{print "Control", "Control", $3,$4,$5,$6}' ${tmpfile}_old.fam > ${tmpfile}_control.fam

  plink --bfile ${tmpfile}_control --bmerge ${tmpfile}_control --merge-mode 7 --out ${tmpfile}_concord

  #get concordance rate from log file
  echo "$ncontrols control samples found" >> $reportfile
  nonmiss=`grep nonmissing ${tmpfile}_concord.log | awk '{print $4}'`  #number non-miss calls
  concord=`grep concordant ${tmpfile}_concord.log | awk '{print $1}'`   #number concordant calls
  rate=`grep concordant ${tmpfile}_concord.log | awk '{print $8}' | cut -d. -f1,2`   #concordance rate (need to remove trailing period)
  echo "Number of non-missing genotype calls:  $nonmiss" >> $reportfile
  echo "Number of concordant calls:  $concord" >> $reportfile
  echo "Concordance rate:  $rate" >> $reportfile
else
  echo "No control samples found" >> $reportfile
fi

rm $tmpfile*
