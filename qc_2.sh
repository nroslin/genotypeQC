#!/bin/bash
module load plink/1.90b3x
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
scriptdir=/hpf/projects/arnold/users/mlemire/scripts/qc/genotypingarrays



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

plink --bfile $dir/$prefix  --extract  ${tmpfile}.common --exclude ${prefix}.exclude.txt  --make-bed --out ${tmpfile}_${prefix} 
#awk '{$5<$6?al=$5":"$6:al=$6":"$5}{$2=$1":"$4; $3=0}{print}' ${tmpfile}_${prefix}.bim | sed 's/ /\t/g'  > ${tmpfile}.tmp
#\mv  ${tmpfile}.tmp  ${tmpfile}_${prefix}.bim 


# pruning LD 
plink --bfile  ${tmpfile}_${prefix} --maf 0.05 --indep-pairwise 1500 100 0.1 --out ${tmpfile}_${prefix}
plink --bfile  ${tmpfile}_${prefix} --extract  ${tmpfile}_${prefix}.prune.in --make-bed --out ${tmpfile}_${prefix}_pruned 

# extracting in 1kg. recoding snps names first 
# awk '{$5<$6?al=$5":"$6:al=$6":"$5}{$2=$1":"$4; $3=0}{print}'  $kgdir/indep.bim > ${tmpfile}_indep.bim 
plink --bed $kgdir/indep.bed --fam $kgdir/indep.fam --bim $kgdir/indep.bim    --extract ${tmpfile}_${prefix}.prune.in  --make-bed --out ${tmpfile}_1kg  

plink --bfile  ${tmpfile}_${prefix} --bmerge ${tmpfile}_1kg.bed ${tmpfile}_1kg.bim  ${tmpfile}_1kg.fam --exclude ${prefix}.exclude.txt --make-bed  --out ${tmpfile}

# first round of flipping 
plink --bed $kgdir/indep.bed --fam $kgdir/indep.fam --bim $kgdir/indep.bim    --extract ${tmpfile}_${prefix}.prune.in --flip  ${tmpfile}-merge.missnp --make-bed --out ${tmpfile}_1kg  
plink --bfile  ${tmpfile}_${prefix}_pruned --bmerge ${tmpfile}_1kg.bed ${tmpfile}_1kg.bim  ${tmpfile}_1kg.fam --exclude ${prefix}.exclude.txt --make-bed  --out ${tmpfile}_merge

plink --bfile ${tmpfile}_merge --pca --out ${prefix}.qc_pca 


R --no-save --args ${prefix}.qc_pca.eigenvec  ${tmpfile}_merge.fam  ${prefix}.1kgpca < $scriptdir/qc_pca.r

\rm ${tmpfile}* 




