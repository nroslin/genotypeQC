#!/bin/bash

module load plink/1.90b6.21
module load R/3.6.1
module load bcftools/1.9


# usage, e.g. ./qc_1.sh ../../original SHA15345
# bfiles are expected to be found there:
dir=$1
# prefix of the bfiles
prefix=$2

# tmp file prefix. All _$$_tmp_* will be deleted
tmpfile=_$$_tmp_ 
# where the scripts called here are located 
scriptdir=/hpf/projects/arnold/users/nroslin/Scripts/QC/GenotypeArrays

imissrate=$3
lmissrate=$4

hetbprange=$5 

echo FID IID SOURCE > ${prefix}.remove.txt 
echo SNP SOURCE > ${prefix}.exclude.txt 

#############
# SAMPLE QC #
#############

##########
# sex check in PLINK sucks, so do it more manually
#use bcftools to count heterozygous calls on X and number of calls Y
#note that PLINK will give warning about non-valid alleles; these are D/I
plink --memory 8000 --bfile $dir/$prefix --chr 23 --recode vcf --snps-only --out _$$_tmp_c23
bcftools +smpl-stats _$$_tmp_c23.vcf > _$$_tmp_c23.bcfout

#change chr Y chromosome code to autosome and calculate missing rate
#otherwise, females will always show no calls
#also, if change sex to missing will produce no results
plink --memory 8000 --bfile $dir/$prefix --chr 24 --make-bed --out _$$_tmp_c24
sed 's/^24/2/' _$$_tmp_c24.bim > _$$_tmp_c24new.bim
plink --memory 8000 --bfile _$$_tmp_c24 --bim _$$_tmp_c24new.bim --missing --out _$$_tmp_c24

#format of bcftools output is awkward to read into R, so format using perl
$scriptdir/qc_sexCheck.pl _$$_tmp $dir $prefix
  #makes $prefix_sexCheck.txt (file with counts and stats, no inference)

#make some plots and do inference
R --no-save --args $prefix < $scriptdir/qc_sexCheck.r > qc_sexCheck.log
  #makes $prefix_inferredSex.txt, $prefix_sexCheck.pdf

#problem samples: inferred sex unknown, OR pedSex not missing and pedSex ne 
#inferred sex
sed '1d' ${prefix}_inferredSex.txt | awk '$4==0 || ( $3!=0 && $3!=$4) { print $1,$2,"SexCheck"}' >> $prefix.remove.txt

# THIS NEEDS TO BE REFINED. AWAITING SEX INFO 
# as of v > 0.2.2 this is commented 
# as os v0.3.0 this is uncommented 
#awk 'NR>1 && $3!=0 && $3!=$4 && $4!=0 && $5=="PROBLEM" {print $1,$2,"SEXCHECK"}' ${prefix}_check-sex.sexcheck  >> ${prefix}.remove.txt 



###########
# creating a new fam file that includes inferred sex when real sex info is not avail 
#also, set any hh genotypes to missing in inferred males

awk '{$3==0?sex=$4:sex=$3} NR>1 {print $1,$2,sex}'    ${prefix}_inferredSex.txt > ${tmpfile}_update_sex 
awk '$3==1 {print $1,$2}' ${tmpfile}_update_sex > ${tmpfile}_males 

plink --memory 8000  --bfile $dir/$prefix --update-sex ${tmpfile}_update_sex  --make-bed --out ${tmpfile}_update_sex --set-hh-missing 




##########
# missingness
plink --memory 8000  --bfile ${tmpfile}_update_sex --missing --out ${prefix}_missing_step1
# excluding SNPs with missing rate > $lmissrate before calculating imiss
awk 'NR>1 && ( $5>'$lmissrate' || $4==0 ) {print $2}' ${prefix}_missing_step1.lmiss > ${tmpfile}.exclude.txt

#find vars in males only on chrX with high missing rate (poorly typed or a lot
#of hh, maybe pseudoautosomal region)
#this is new as of v0.2.4 
plink --memory 8000  --bfile ${tmpfile}_update_sex --missing --chr 23 --out ${prefix}_missing_Xmales --keep ${tmpfile}_males
awk 'NR>1 && ( $5>'$lmissrate' || $4==0 ) {print $2}' ${prefix}_missing_Xmales.lmiss >> ${tmpfile}.exclude.txt


# this file not needed
\rm ${prefix}_missing_step1.imiss
#remove these X chrom vars before calculating missing rate per sample
plink --memory 8000  --bfile ${tmpfile}_update_sex --missing --exclude ${tmpfile}.exclude.txt --out ${prefix}_missing_step2
# this file  not needed 
\rm ${prefix}_missing_step2.lmiss
awk 'NR>1 && $6>'$imissrate' {print $1,$2,"IMISS"}' ${prefix}_missing_step2.imiss  >> ${prefix}.remove.txt 


##########
# het
awk '$1<23 {print $2}' $dir/${prefix}.bim > ${tmpfile}.extract.txt
plink --memory 8000  --bfile ${tmpfile}_update_sex --het --extract ${tmpfile}.extract.txt --exclude  ${tmpfile}.exclude.txt --out ${prefix}_het
# this excludes samples if F is and outlier wrt boxplot with range=$hetbprange
Rscript $scriptdir/qc_het.r ${prefix}_het.het ${prefix}.remove.txt $hetbprange

#############
# MARKER QC #
#############


##############
# missingness

plink --memory 8000  --bfile ${tmpfile}_update_sex --remove ${prefix}.remove.txt --missing --out ${prefix}_missing_step3
awk 'NR>1 && ( $5>'$lmissrate' ) {print $2,"LMISS"}' ${prefix}_missing_step3.lmiss >> ${prefix}.exclude.txt 

# chrX-specific
plink --memory 8000  --bfile ${tmpfile}_update_sex --missing --chr 23 --out ${prefix}_missing_Xmales_step3 --keep ${tmpfile}_males --remove  ${prefix}.remove.txt 
awk 'NR>1 && ( $5>'$lmissrate' ) {print $2,"LMISS"}' ${prefix}_missing_Xmales_step3.lmiss >> ${prefix}.exclude.txt

\rm ${prefix}_missing_step3.imiss



#################
# duplicated SNPs
# ignoring SNPs with 0 allele, as both SNPs may refer to different alleles 

# recoding SNP names with chr:pos:a1:a2 ; alleles are ordered 
# IGNORING THE ALLELE SINCE THE STRAND MAY BE DIFFERENT -- JUST FOCUSING ON POSITION
# AS OF v0.2.4 I am not including SNPs with chr==0 among duplicates 

#awk '{$5<$6?al=$5":"$6:al=$6":"$5}{print $1":"$4":"al, $2}' ${tmpfile}_update_sex.bim | sort -k 1,1 > ${tmpfile}.snpnames  
awk '{$5<$6?al=$5":"$6:al=$6":"$5} $1>0 {print $1":"$4, $2}' ${tmpfile}_update_sex.bim | sort -k 1,1 > ${tmpfile}.snpnames  


# extracting SNP names with same pos and alleles 

#cut -f 1 -d' '  ${tmpfile}.snpnames  | sort | uniq -c | awk '$1>1{print $2}' |\
#  awk '{split( $1,pos,":")} pos[3]!=0 && pos[4]!=0 {print}'  | sort -k 1,1 > ${tmpfile}.duplicates 

cut -f 1 -d' '  ${tmpfile}.snpnames  | sort | uniq -c | awk '$1>1{print $2}' |\
   sort -k 1,1 > ${tmpfile}.duplicates 


join -1 1 -2 1 ${tmpfile}.duplicates  ${tmpfile}.snpnames | sort -k 2,2 >  ${tmpfile}.snpnames.sorted 

# want to be picking the one with most call rate 
awk '{print $2,$5}'  ${prefix}_missing_step3.lmiss | sort -k 1,1 > ${tmpfile}.lmiss 

# the sort -r is so that I will tend to pick the rs designation if ties but won't always work 
join -1 2 -2 1 ${tmpfile}.snpnames.sorted   ${tmpfile}.lmiss |\
 awk '{print $2,$1,$3}' | sort -r |\
 awk 'BEGIN{last=""; smallest=999; best="" }{
   if($1!=last){ if( NR>1){print best}; smallest=$3; best=$2; last=$1 }
   else { if( $3< smallest ){ best=$2; smallest=$3 } } }' | sort  > ${tmpfile}.bestcallrate 


# bug fixed in v0.2.1: used to print $2
join -v 2 -1 1 -2 2  ${tmpfile}.bestcallrate   ${tmpfile}.snpnames.sorted   | awk '{print $1,"DUPLICATE"}'  >>  ${prefix}.exclude.txt 

sort -u  ${prefix}.exclude.txt > ${tmpfile} ; mv ${tmpfile} ${prefix}.exclude.txt



########
# chr0
awk '$1==0 {print $2,"CHR_0"}' $dir/$prefix.bim >> ${prefix}.exclude.txt






\rm ${tmpfile}* 



















#############################################################################
#   Copyright 2019 Mathieu Lemire
#             2022 Nicole M. Roslin
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
