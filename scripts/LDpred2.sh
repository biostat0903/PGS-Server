#!/bin/bash

while getopts ":s:G:P:c:C:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    G) val_geno="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    C) chr="$OPTARG"
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
if [ "${cov}" -eq 0 ]
then 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"

CURRENT_DIR=$(cd `dirname $0`; pwd)
PACKAGE_DIR=/home/yasheng/comprsWeb/scripts/
LDPRED2CODE=${PACKAGE_DIR}/LDpred2.R

if [ "${cov}" -eq 0 ]
then 
echo 0000
Rscript ${LDPRED2CODE} --summ ${summary_file_prefix} --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                       --dist 1000 --chr ${chr} --outpath ${outpath}
else
Rscript ${LDPRED2CODE} --summ ${summary_file_prefix} --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                       --dist 1000 --chr ${chr} --covariates ${cov} --outpath ${outpath}
fi 