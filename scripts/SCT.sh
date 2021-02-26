#!/bin/bash

while getopts ":s:G:P:p:r:d:c:t:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    G) val_geno="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    p) p_len="$OPTARG"
    ;;
    r) r2_val="$OPTARG"
    ;;
    d) dist_str="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    t) thread="$OPTARG"
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
printf "\033[33mArgument p_len is %s  \033[0m\n" "$p_len"
printf "\033[33mArgument r2_val is %s  \033[0m\n" "$r2_val"
printf "\033[33mArgument dist_str is %s  \033[0m\n" "$dist_str"
if [ ${cov} -ne 0 ]
then 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"


CURRENT_DIR=$(cd `dirname $0`; pwd)
PACKAGE_DIR=/home/yasheng/comprsWeb/scripts/
SCTCODE=${PACKAGE_DIR}/SCT.R

if [ ${cov} -eq 0 ]
then 
Rscript ${SCTCODE} --summ ${summary_file_prefix} --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                   --p_len ${p_len} --r2_val ${r2_val} --dist_str ${dist_str} --thread ${thread} --outpath ${outpath}
else 
Rscript ${SCTCODE} --summ ${summary_file_prefix} --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                   --p_len ${p_len} --r2_val ${r2_val} --dist_str ${dist_str} --covariates ${cov} --thread ${thread} --outpath ${outpath}
fi
