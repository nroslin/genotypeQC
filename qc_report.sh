#!/bin/bash

#############################
# args

dir=$1 
prefix=$2


###################################
# setting and grepping variables 

scriptdir=/hpf/projects/arnold/users/mlemire/scripts/qc/genotypingarrays

today=`date`
here=`pwd -P`

# only after commit f254c768e01f that these are set in README.steps 
imissrate=`grep imissrate README.steps | head -1 | cut -f 2 -d'='`
lmissrate=`grep imissrate README.steps | head -1 | cut -f 2 -d'='`

# legacy, default was 0.03 hard coded in qc_1.sh 
if [ -z $imissrate ]; then
  imissrate=0.03 
  echo \!\!\!\!\!\!\!\!\! CONFIRM imissrate=$imissrate
fi
if [ -z $lmissrate ]; then
  lmissrate=0.03 
  echo \!\!\!\!\!\!\!\!\! CONFIRM lmissrate=$lmissrate
fi





###############################
# writing report 

echo "0        1         2         3         4         5         6         7"
echo 1234567890123456789012345678901234567890123456789012345678901234567890
echo "######################################################################" 
echo "#" Report produced on $today 
echo "#" $here
echo "######################################################################" 
echo 
echo "   ---------------------"
echo "1. Sample and SNP counts"
echo "   ---------------------"

nsamples=`wc -l $dir/$prefix.fam | awk '{print $1}'`
ncontrols=`cat *.remove.txt | awk '$3=="CONTROL" {print}' | sort -u |wc -l | awk '{print $1}'`
nsnps=`wc -l $dir/$prefix.bim | awk '{print $1}'`
ndups=`cat *.exclude.txt | awk '$2=="DUPLICATE" {print}' | sort -u |wc -l | awk '{print $1}'`

echo "Number of DNA samples genotyped: $nsamples"
if [ $ncontrols -gt 0 ]; then 
 echo "   including $ncontrols internal controls"
 echo "   listed in the following file(s):"
 files=`grep CONTROL *.remove.txt | cut -f 1 -d':' | sort -u`
 for f in $files; do 
  echo "      " $here/$f
 done 
fi 
echo 
echo "Number of SNPs genotyped: $nsnps"
if [ $ndups -gt 0 ]; then 
 echo "   (including $ndups duplicated SNPs)"
fi 



echo 
echo "   --------------"
echo "2. SNP statistics" 
echo "   --------------"

# average snp call rate 
snpavgcr=`awk 'BEGIN{tot=0;n=0} NR>1 {tot=tot+$5; n=n+1} END {print (1-tot/n)}' ${prefix}_missing_step1.lmiss`
snpfailcr=`grep LMISS ${prefix}.exclude.txt |wc -l |awk '{print $1}'`


medianfailedcr=`awk -v m=$lmissrate 'NR>1 && $5>m  {print 1-$5}' ${prefix}_missing_step1.lmiss |\
 sort -g | awk -v n=$snpfailcr 'NR>=n/2 {print int( 1000*$1+.5 )/10}' | head -1`



echo 
echo "Average SNP call rate: $snpavgcr"

echo 
echo "Number of SNPs with call rate less than "`echo $lmissrate | awk '{print (1-$1)}'`": $snpfailcr"
if [ $snpfailcr -gt 0 ]; then 
 medianfailedcr=`awk -v m=$lmissrate 'NR>1 && $5>m  {print 1-$5}' ${prefix}_missing_step1.lmiss |\
                     sort -g | awk -v n=$snpfailcr 'NR>=n/2 {print $1}' | head -1`
 echo " median call rate among them: ${medianfailedcr}"
fi

# Hardy

noutliers=`wc -l ${prefix}.1kgpca.outliers.txt |awk '{print $1}'`
nhardy=`grep HWE: ${prefix}.exclude.txt | wc -l | awk '{print $1}'`

echo 
echo "Number of SNPs that failed HWE at FDR<0.01: $nhardy"
echo " (calculated after removal of $noutliers outlying samples based on ancestry)" 





echo 
echo "   -----------------"
echo "3. Sample statistics" 
echo "   -----------------"

# average sample call rate 
sampleavgcr=`awk 'BEGIN{tot=0;n=0} NR>1 {tot=tot+$6; n=n+1} END {print (1-tot/n)}' ${prefix}_missing_step2.imiss`
samplefailcr=`grep IMISS ${prefix}.remove.txt |wc -l |awk '{print $1}'`


echo "Average sample call rate: $sampleavgcr"
echo " (among SNPs with call rate above "`echo $lmissrate | awk '{print (1-$1)}'`")"
echo 
echo "Number of samples with call rate less than "`echo $imissrate | awk '{print (1-$1)}'`": $samplefailcr"
echo " (among SNPs with call rate above "`echo $lmissrate | awk '{print (1-$1)}'`")"
if [ $samplefailcr -gt 0 ]; then 
 medianfailedcr=`awk -v m=$imissrate 'NR>1 && $6>m  {print 1-$6}' ${prefix}_missing_step2.imiss |\
                     sort -g | awk -v n=$samplefailcr 'NR>=n/2 {print $1}' | head -1`
 echo " median call rate among them: ${medianfailedcr}"
fi



wrongsex=`grep PROBLEM ${prefix}_check-sex.sexcheck | awk '$3!=0 {print}' |wc -l | awk '{print $1}'`
echo
echo "Number of samples with wrong sex: $wrongsex"


nhet=`awk '$3=="HET" {print}' ${prefix}.remove.txt |wc -l | awk '{print $1}'`
echo 
echo "Number of samples with excess/deficit of heterozygosity: $nhet"


echo
echo "   ---------------------------"
echo "4. Technical exclusion summary"
echo "   ---------------------------"


nsnpsexcl=`expr $snpfailcr + $nhardy`

echo
echo "Number of SNPs excluded (ignoring duplicates): $nsnpsexcl / $nsnps ("`echo $nsnpsexcl $nsnps | awk '{print int(1000*$1/$2+.5)/10}'`"%)"


nsamplesexcl=`expr $samplefailcr + $wrongsex +  $nhet `
echo 
echo "Number of samples excluded: $nsamplesexcl / $nsamples ("`echo $nsamplesexcl $nsamples | awk '{print int(1000*$1/$2+.5)/10}'`"%)"

echo
echo 

















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
