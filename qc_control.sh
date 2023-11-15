#!/bin/bash

# Created 15 November 2023.
# Last modified:  15 Nov 2023

# Process control samples.  Assume if ID does not start with TAG then is a
# control.  Also assume all controls are the same ID (NA12878) and so should
# have identical genotypes.

module load plink/1.90b6.21


echo "#####################################"
echo "Calling script $0"
echo "Processing of control samples"

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
echo "Processing of control samples" >> $reportfile

awk '$2!~/^TAG/ {print $1,$2}' $dir/$prefix.fam > $tmpfile.control
ncontrols=`wc -l $tmpfile.control | awk '{print NF}'`

if [ $ncontrols -gt 0 ]
then
  #there are controls
  plink --bfile $dir/$prefix --exclude $prefix.exclude.txt --keep $tmpfile.control --genome --out $tmpfile

  #expecting all pairs to have hat(pi) = 1
  echo "$ncontrols control samples found" >> $reportfile
  sed '1d' $tmpfile.genome | awk '$10!=1 {print $1,$2,$3,$4,$10}' > $tmpfile.nodup
  echo "The following control pairs do not have duplicate genotypes" >> $reportfile
  nlines=`wc -l $tmpfile.nodup | awk '{print $1}'`
  if [ $nlines -gt 0 ]
  then
	echo "FID1 IID1 FID2 IID2 PI_HAT" >> $reportfile
	cat $tmpfile.nodup
  else
	echo "NONE" >> $reportfile
  fi

else
  echo "No control samples found" >> $reportfile
fi

rm $tmpfile*
