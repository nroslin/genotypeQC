#!/bin/bash
#SBATCH --tmp=8g

module load plink/1.90b6.21
module load R/3.5.1

#alias plink='plink --noweb'

# 1kg plink files provided by TCAG
# avail at http://www.tcag.ca/tools/1000genomes.html 
kgdir=/hpf/projects/arnold/references/www.tcag.ca/documents/tools




# bfiles are expected to be found there:
dir=$1
# prefix of the bfiles
prefix=$2
#should we include AMR in ancestry inference?  1=yes, 0=no
amr=$3
#do we want to detect outliers in the samples?  1=yes, 0=no
outlier=$4

reportfile=${prefix}_QCreport.txt   #file was started in qc_1.sh

# tmp file prefix. All _$$_tmp_* will be deleted
tmpfile=_$$_tmp_ 
# where the scripts called here are located 
scriptdir=/hpf/projects/arnold/users/nroslin/Scripts/Genotypes/qc

echo "#################################################################"
echo "Running script $scriptdir/$0"
echo "Command-line arguments:"
echo "PLINK binary files directory:  $1"
echo "PLINK files prefix:  $2"
echo "Ancestry inference including AMR?  $3"
echo "Ancestry outlier detection?  $4"
echo
echo


awk '{print $2,$5,$6}' $kgdir/indep.bim  | sort -k 1,1 > ${tmpfile}.1kg
awk '$1<=22 {print $2,$5,$6}' $dir/${prefix}.bim  | sort -k 1,1 > ${tmpfile}.bim

# excluding A/T and C/G SNPs 
# common based on chr:pos:a1:a2 
join -1 1 -2 1 ${tmpfile}.bim ${tmpfile}.1kg |\
 awk 'BEGIN{s["A"]="T";s["C"]="G";s["G"]="C";s["T"]="A" }
      (($2==$4 && $3==$5 ) || ($3==$4 && $2==$5 ) || ($2==s[$4] && $3==s[$5] ) || ($3==s[$4] && $2==s[$5] )) && 
      $2$3!="AT" && $2$3!="TA" && $2$3!="CG" && $2$3 !="GC" {print $1}'  > ${tmpfile}.common 


# extracting first in one dataset because plink complains about strand for SNPs not in the extract list 
# also recoding the SNP name based on chr:pos
# duplicates should have been flagged in the exclude file 

plink --memory 8000 --bfile $dir/$prefix  --extract  ${tmpfile}.common --exclude ${prefix}.exclude.txt  --make-bed --out ${tmpfile}_${prefix} 
#awk '{$5<$6?al=$5":"$6:al=$6":"$5}{$2=$1":"$4; $3=0}{print}' ${tmpfile}_${prefix}.bim | sed 's/ /\t/g'  > ${tmpfile}.tmp
#\mv  ${tmpfile}.tmp  ${tmpfile}_${prefix}.bim 


# pruning LD 
plink --memory 8000 --bfile  ${tmpfile}_${prefix} --maf 0.05 --indep-pairwise 1500 100 0.2 --out ${tmpfile}_${prefix}
plink --memory 8000 --bfile  ${tmpfile}_${prefix} --extract  ${tmpfile}_${prefix}.prune.in --make-bed --out ${tmpfile}_${prefix}_pruned 


############## 
# looking for duplicates
# want to do this before PCA in case there are twins
# eg. the duplicated control samples will affect the PCA results, and twin
# pairs tend to be flagged as outliers

# This section looks for samples who are expected to be duplicates (IIDs are the same)
# Removed if they do not have identical genotypes (FALSE_DUPLICATES)
# Otherwise, keep the one with higher call rate

plink --memory 8000 --bfile ${tmpfile}_${prefix}_pruned --genome full --out ${tmpfile}_${prefix}_pruned

awk '$2==$4 && $10< 0.9 {print $1,$2,"FALSE_DUPLICATES\n"$3,$4,"FALSE_DUPLICATES"}' ${tmpfile}_${prefix}_pruned.genome >> ${prefix}.remove.txt 

awk '$2==$4 && $10 > 0.9 {print $1,$2"\n"$3,$4}' ${tmpfile}_${prefix}_pruned.genome > ${prefix}.duplicates.txt 

#if there are duplicates, calculate call rate and keep the one with the
#highest CR
ndups=`wc -l ${prefix}.duplicates.txt | awk '{print $1}'`
if [ $ndups -gt 0 ]
then
  plink --bfile $dir/$prefix --missing --out ${prefix}.duplicates --keep ${prefix}.duplicates.txt
  dups=`awk '{print $2}' ${prefix}.duplicates.txt | sort -u`
   for d in $dups; do 
	awk '$2=="'$d'" {print}' ${prefix}.duplicates.imiss | sort -k 6,6 -g | tail -n +2 | awk '{print $1,$2,"DUPLICATE_HIGHER_MISS"}' >> ${prefix}.remove.txt
  done 
fi

### report ###
echo >> $reportfile
echo >> $reportfile
echo "Sample duplicates" >> $reportfile
echo "-----------------" >> $reportfile
echo "The following samples were removed because they were expected to be part of a duplicate/replicate but are not:" >> $reportfile
lines=`grep FALSE_DUPLICATES ${prefix}.remove.txt | wc -l | awk '{print $1}'`
if [ $lines -gt 0 ]
then
  grep FALSE_DUPLICATES ${prefix}.remove.txt | awk '{print $1,$2}' >> $reportfile
else
  echo "NONE" >> $reportfile
fi
echo >> $reportfile

echo "The following samples were removed because they were part of a known duplicate/replicate and had lower call rate:" >> $reportfile
lines=`grep DUPLICATE_HIGHER_MISS ${prefix}.remove.txt | wc -l | awk '{print $1}'`
if [ $lines -gt 0 ]
then
  grep DUPLICATE_HIGHER_MISS ${prefix}.remove.txt | awk '{print $1,$2}' >> $reportfile
else
  echo "NONE" >> $reportfile
fi
echo >> $reportfile

echo "The following sample duplicates remain.  They could be known twins or" >> $reportfile
echo "accidental duplicates:  check expected twin status" >> $reportfile
sed '1d' ${tmpfile}_${prefix}_pruned.genome | awk '$10>0.9 && $2!=$4 {print $1,$2}' > ${tmpfile}.mz   #1st ID in pair
sed '1d' ${tmpfile}_${prefix}_pruned.genome | awk '$10>0.9 && $2!=$4 {print $3,$4}' >> ${tmpfile}.mz   #2nd ID in pair
lines=`wc -l ${tmpfile}.mz | awk '{print $1}'`
if [ $lines -gt 0 ]
then
  echo "FID IID" >> $reportfile
  sort -u ${tmpfile}.mz >> $reportfile
else
  echo "NONE" >> $reportfile
fi
echo >> $reportfile
### report ###

#find clusters of duplicates, and assign them to sibships of appropriate size
#using code from Mathieu Lemire
#we don't care if this is expected or not, we just want to remove dups before
#doing PCA
sed '1d' ${tmpfile}_${prefix}_pruned.genome | awk '$10>0.9 {print $1 "_" $2, $3 "_" $4}' | sort -k1,1 > _mz1   #find all pairs
#duplicate each pair, with different person listed first
awk '{print $1,$2"\n"$2,$1}' _mz1 | sort -k1,1 > _mzpairs

#building up clusters of duplicates
join -1 1 -2 1 _mzpairs  _mzpairs  |\
 awk '{printline=1}{for(i=1; i<= NF-1; i++ ){if($i>=$(i+1)){printline=0} }}
      printline==1 {print}'  | sort -k 1,1 > _tmp_mztwins

continue=`wc -l _tmp_mztwins | awk '{print $1}'`
i=1
while [ $continue -gt 0 ]; do
 i=$[$i+1]
 mv _tmp_mztwins _mz$i
 join -1 1 -2 1 _mz$i _mzpairs |\
  awk '{printline=1}{for(i=1; i<= NF-1; i++ ){if($i>=$(i+1)){printline=0} }}
       printline==1 {print}'  | sort -k 1,1 > _tmp_mztwins
 continue=`wc -l _tmp_mztwins | awk '{print $1}'`
done

# Fx_n where x is mzhip size, n is current number of mzhip of size x 
# Fx_n where x is mzhip size, n is current number of mzhip of size x 

rm -f  _mztwins ; touch  _mztwins
while [ $i -gt 0 ]; do
 awk '{for(j=1; j<=NF-1; j++ ){print $j}}' _mztwins | sort -u > _mz_remove
 awk '{for(j=1; j<=NF; j++ ){print $j,"MZ""'$i'"+1"_"NR}}' _mz$i | sort > _tmp_mztwins
 join -v 2 -1 1 -2 1 _mz_remove  _tmp_mztwins  >> _mztwins
 i=$[$i-1]
done

# check that this is empty; all samples are assigned a unique mztwins 
echo 
echo "Twin check:  there should be nothing immediately following this line"
sort -u _mztwins | awk '{print $1}'  |sort | uniq -c | awk '$1>1{print}'
echo
echo

sort -k 1,1  _mztwins > ${prefix}_mztwins.txt
\rm _mz* _tmp_mz*

#for PCA, select one duplicate to keep and remove the rest
#find unique cluster IDs
awk '{print $2}' ${prefix}_mztwins.txt | sort -u > ${tmpfile}.clustid
#remove everyone in cluster except first ID
rm -rf ${tmpfile}.toRemoveTwins.txt; touch ${tmpfile}.toRemoveTwins.txt
for clustid in `cat ${tmpfile}.clustid`
do
  grep -w $clustid ${prefix}_mztwins.txt | sed '1d' | awk '{print $1}' | sed 's/_/ /' >> ${tmpfile}.toRemoveTwins.txt
done
nsamples=`wc -l ${tmpfile}.toRemoveTwins.txt | awk '{print $1}'`
echo "The following $nsamples duplicates will be removed from PCA analysis."
echo "They will be reported in ${prefix}.1kgpca.closestAncestry.txt and ${prefix}.1kgpca.outliers.txt with the same information as their twin."
echo
cat ${tmpfile}.toRemoveTwins.txt
echo
echo

#note that this file is just to perform a reasonable PCA, and is not the final
#list of IDs to remove
#remove these IDs from observed data
plink --memory 8000 --bfile  ${tmpfile}_${prefix}_pruned --remove ${tmpfile}.toRemoveTwins.txt --make-bed  --out ${tmpfile}_${prefix}_pruned2

############## 
# Add in 1kg samples and do PCA

# extracting in 1kg. recoding snps names first 
# awk '{$5<$6?al=$5":"$6:al=$6":"$5}{$2=$1":"$4; $3=0}{print}'  $kgdir/indep.bim > ${tmpfile}_indep.bim 

### AMR
if [ $amr -eq 0 ]
then
  ### no AMR:  exclude any AMR samples
  #note that FID is not necessarily correct, so get from fam file
  awk '$4=="AMR" {print $2}' $kgdir/sampleTable.txt | egrep -f - $kgdir/indep.fam | awk '{print $1,$2}' > ${tmpfile}_amr.txt
  #also remove 1kg population outliers (listed in 1kg report)
  cat ${tmpfile}_amr.txt $kgdir/popOutliers.txt | sort -u > ${tmpfile}_1kg.remove.txt
else
  cp $kgdir/popOutliers.txt ${tmpfile}_1kg.remove.txt
fi

plink --memory 8000 --bed $kgdir/indep.bed --fam $kgdir/indep.fam --bim $kgdir/indep.bim --remove ${tmpfile}_1kg.remove.txt   --extract ${tmpfile}_${prefix}.prune.in  --make-bed --out ${tmpfile}_1kg  

plink --memory 8000 --bfile  ${tmpfile}_${prefix} --bmerge ${tmpfile}_1kg.bed ${tmpfile}_1kg.bim  ${tmpfile}_1kg.fam --exclude ${prefix}.exclude.txt --make-bed  --out ${tmpfile}

# first round of flipping 
plink --memory 8000 --bed $kgdir/indep.bed --fam $kgdir/indep.fam --bim $kgdir/indep.bim  --remove ${tmpfile}_1kg.remove.txt  --extract ${tmpfile}_${prefix}.prune.in --flip  ${tmpfile}-merge.missnp --make-bed --out ${tmpfile}_1kg  
plink --memory 8000 --bfile  ${tmpfile}_${prefix}_pruned2 --bmerge ${tmpfile}_1kg.bed ${tmpfile}_1kg.bim  ${tmpfile}_1kg.fam --exclude ${prefix}.exclude.txt --make-bed  --out ${tmpfile}_merge --allow-no-sex

plink --memory 8000 --bfile ${tmpfile}_merge --pca --out ${prefix}.qc_pca 


mv ${tmpfile}_merge.fam  ${prefix}.qc_pca.fam
rm -f ${prefix}.1kgpca.outliers.txt
R --no-save --args ${prefix}.qc_pca.eigenvec  ${prefix}.qc_pca.fam  ${prefix}.1kgpca $amr $outlier < $scriptdir/qc_pca.r




#need to go back and restore information of twins who were removed (assign
#them the same ancestry as their twin)

lines=`wc -l ${tmpfile}.toRemoveTwins.txt | awk '{print $1}'`
if [ $lines -gt 0 ]
then
  #add information for removed twins into closestAncestry file
  awk '{print $2}' ${tmpfile}.toRemoveTwins.txt > ${tmpfile}.rmid  #removed IDs
  for id in `cat ${tmpfile}.rmid`
  do
	twinid=`grep $id ${prefix}_mztwins.txt | awk '{print $2}'`  #id of this twin
	withdata=`grep $twinid ${prefix}_mztwins.txt | head -1 | awk '{print $1}' | cut -f2 -d_`   #the first listed twin should have data
	grep $withdata ${prefix}.1kgpca.closestAncestry.txt | cut -f3- -d" " > ${tmpfile}.mydata   #data for the analyzed twin, stripping off FID, IID
	grep $id ${prefix}_mztwins.txt | awk '{print $1}' | sed 's/_/ /' > ${tmpfile}.myid
	paste -d" " ${tmpfile}.myid ${tmpfile}.mydata >> ${prefix}.1kgpca.closestAncestry.txt

	if [ -e ${prefix}.1kgpca.outliers.txt ]
	then
	  outlier=`grep $withdata ${prefix}.1kgpca.outliers.txt | wc -l`
	  if [ $outlier -gt 0 ]
	  then
		echo "Sample $withdata is an outlier, twin $id added to outlier list"
		cat ${tmpfile}.myid >> ${prefix}.1kgpca.outliers.txt
		echo
	  fi
	fi
  done

  ##note:  need to check if the twin is in the outlier list as well
fi



\rm ${tmpfile}* 













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
