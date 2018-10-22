#!/bin/bash
module load plink/1.07
module load R/3.5.1

alias plink='plink --noweb'

# 1kg plink files provided by TCAG
kgdir=/hpf/projects/arnold/references/www.tcag.ca/documents/tools

# bfiles are expected to be found there:
dir=$1
# prefix of the bfiles
prefix=$2

# tmp file prefix. All _$$_tmp_* will be deleted
tmpfile=_$$_tmp_ 
# where the scripts called here are located 
scriptdir=/hpf/projects/arnold/data/genotypes/Redo_-SHA15345_PSYCHIP_PLATE_1-5/qc/v0.1


awk '{print $1":"$4,$5,$6}' $kgdir/snpsForPCA.bim  | sort -k 1,1 > ${tmpfile}.1kg
awk '{print $1":"$4,$2,$5,$6}' $dir/${prefix}.bim  | sort -k 1,1 > ${tmpfile}.bim

# excluding A/T and C/G SNPs 
# common based on chr:pos:a1:a2
join -1 1 -2 1 ${tmpfile}.bim ${tmpfile}.1kg |\
 awk '( ($3==$5 && $4==$6 ) || ($4==$5 && $3==$6 ) ) && $3$4!="AT" && $3$4!="TA" && $3$4!="CG" && $3$4 !="GC" {print $2}'  > ${tmpfile}.common 

# merging
# extracting first in one dataset because plink complains about strand for SNPs not in the extract list 
# also recoding the SNP name 
plink --bfile $dir/$prefix  --extract  ${tmpfile}.common --exclude ${prefix}.exclude.txt  --make-bed --out ${tmpfile}_${prefix} 
awk '{$5<$6?al=$5":"$6:al=$6":"$5}{$2=$1":"$4":"al; $3=0}{print}' ${tmpfile}_${prefix}.bim | sed 's/ /\t/g'  > ${tmpfile}.tmp
\mv  ${tmpfile}.tmp  ${tmpfile}_${prefix}.bim 

# need to recode snp name based on position
awk '{print $2}' ${tmpfile}_${prefix}.bim > ${tmpfile}.common 
awk '{$5<$6?al=$5":"$6:al=$6":"$5}{$2=$1":"$4":"al; $3=0}{print}'  $kgdir/indep.bim > ${tmpfile}_indep.bim 
plink --bed $kgdir/indep.bed --fam $kgdir/indep.fam --bim  ${tmpfile}_indep.bim    --extract  ${tmpfile}.common --make-bed --out ${tmpfile}_1kg  

plink --bfile  ${tmpfile}_${prefix} --bmerge ${tmpfile}_1kg.bed ${tmpfile}_1kg.bim  ${tmpfile}_1kg.fam --exclude ${prefix}.exclude.txt --recodeA --out ${tmpfile} 

R --no-save --args ${tmpfile}.raw ${prefix}.1kgpca < $scriptdir/qc_pca.r

\rm ${tmpfile}* 




