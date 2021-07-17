#!/bin/bash

while getopts ":C:s:G:P:p:c:o:t:" opt; do
  case $opt in
    C) codepath="$OPTARG"
    ;;
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
    t) thread="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument codepath is %s  \033[0m\n" "$codepath"
printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument valid_geno is %s  \033[0m\n" "$val_geno"
printf "\033[33mArgument valid_pheno is %s  \033[0m\n" "$val_pheno"
printf "\033[33mArgument pop is %s  \033[0m\n" "$pop"
if [ ${cov} ]
then 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
if [ ${thread} ]
then 
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
else
thread=1
fi
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"


LASSOSUM=${codepath}/lassosum.R

## sample size of GWAS
nobs=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)

if [ !${cov}  ]
then 

Rscript ${LASSOSUM} --summ ${summary_file_prefix}.assoc.txt --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                    --n ${n} --population ${pop} --outpath ${outpath} --thread ${thread}
else

Rscript ${LASSOSUM} --summ ${summary_file_prefix}.assoc.txt --valid_genotype ${val_geno} --valid_phenotype ${val_pheno}\
                    --n ${n} --population ${pop} --covariates ${cov} --outpath ${outpath} --thread ${thread}
fi
