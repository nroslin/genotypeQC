#!/bin/bash

module load plink/1.90b3x
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


plink --memory 8000  --bfile $dir/$prefix --remove ${tmpfile}.remove --hardy --out ${prefix}_hardy 

R --no-save --args ${prefix}_hardy.hwe ${prefix}.exclude.txt  < ${scriptdir}/qc_hwe.r 


\rm ${tmpfile}* 
\rm ${prefix}*.nosex
\rm ${prefix}*.nof
\rm ${prefix}*.hh






















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
