#!/bin/bash

while getopts ":s:G:P:p:c:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    G) val_geno="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument valid_geno is %s  \033[0m\n" "$val_geno"
printf "\033[33mArgument valid_pheno is %s  \033[0m\n" "$val_pheno"
if [ ${cov} -ne 0 ]
then 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"


CURRENT_DIR=$(cd `dirname $0`; pwd)
PACKAGE_DIR=/home/yasheng/comprsWeb/scripts/
LASSOSUM=${PACKAGE_DIR}/lassosum.R

## sample size of GWAS
nobs=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)

if [ ${cov} -eq 0 ]
then 

Rscript ${LASSOSUM} --summ ${summary_file_prefix}.assoc.txt --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                    --n ${n} --population ${pop} --outpath ${outpath}
else

Rscript ${LASSOSUM} --summ ${summary_file_prefix}.assoc.txt --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                    --n ${n} --population ${pop} --covariates ${cov} --outpath ${outpath}
fi