#!/bin/bash
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

# tmp file prefix. All _$$_tmp_* will be deleted
tmpfile=_$$_tmp_ 
# where the scripts called here are located 
scriptdir=/hpf/projects/arnold/users/nroslin/Scripts/Genotypes/qc



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
plink --memory 8000 --bfile  ${tmpfile}_${prefix} --maf 0.05 --indep-pairwise 1500 100 0.1 --out ${tmpfile}_${prefix}
plink --memory 8000 --bfile  ${tmpfile}_${prefix} --extract  ${tmpfile}_${prefix}.prune.in --make-bed --out ${tmpfile}_${prefix}_pruned 

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
plink --memory 8000 --bfile  ${tmpfile}_${prefix}_pruned --bmerge ${tmpfile}_1kg.bed ${tmpfile}_1kg.bim  ${tmpfile}_1kg.fam --exclude ${prefix}.exclude.txt --make-bed  --out ${tmpfile}_merge --allow-no-sex

plink --memory 8000 --bfile ${tmpfile}_merge --pca --out ${prefix}.qc_pca 


mv ${tmpfile}_merge.fam  ${prefix}.qc_pca.fam
R --no-save --args ${prefix}.qc_pca.eigenvec  ${prefix}.qc_pca.fam  ${prefix}.1kgpca $amr < $scriptdir/qc_pca.r


############## 
# duplicated samples 

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
