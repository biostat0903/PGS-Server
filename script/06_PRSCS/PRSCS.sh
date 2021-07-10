#!/bin/bash

while getopts ":C:s:L:T:G:P:o:i:c:" opt; do
  case $opt in
    C) codepath="$OPTARG"
    ;;
    s) summary_file_prefix="$OPTARG"
    ;;
    L) LDpath="$OPTARG"
    ;;
    T) type="$OPTARG"
    ;;
    G) val_geno="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    i) index="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument codepath is %s  \033[0m\n" "$codepath"
printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument type is %s  \033[0m\n" "$type"
printf "\033[33mArgument val_geno is %s  \033[0m\n" "$val_geno"
if [[ "${type}" == "tuning" ]]
then
printf "\033[33mArgument val_pheno is %s  \033[0m\n" "$val_pheno"
printf "\033[33mArgument index is %s  \033[0m\n" "$index"
fi
if [ -n "$cov" ]
then
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"

PRSCS=${codepath}PRScs/PRScs.py
PLINK=${codepath}/plink

awk '{print $2,$6,$7,$9,$11}' ${summary_file_prefix}.assoc.txt > ${summary_file_prefix}.prscs
sed -i '1i\SNP A1 A2 BETA P' ${summary_file_prefix}.prscs
nobs=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)


if [ $type == "auto" ]
then 
for chr in `seq 1 22`
do
python3 ${PRSCS} --ref_dir=${LDpath} --bim_prefix=${val_geno} \
                 --sst_file=${summary_file_prefix}.prscs --n_gwas=${n} \
                 --chrom=${chr} --out_dir=${outpath}/summary --n_iter=10 --n_burnin=5
mv ${outpath}/summary_pst_eff_a1_b0.5_phiauto_chr${chr}.txt ${outpath}/esteff_auto_chr${chr}.txt
done
else 
for chr in `seq 1 22`
do
for phi in 1e-06 1e-04 1e-02 1e+00
do
python3 ${PRSCS} --ref_dir=${LDpath} --bim_prefix=${val_geno} \
                 --sst_file=${summary_file_prefix}.prscs --n_gwas=${n} --phi ${phi}\
                 --chrom=${chr} --out_dir=${outpath}/summary --n_iter=10 --n_burnin=5
mv ${outpath}/summary_pst_eff_a1_b0.5_phi${phi}_chr${chr}.txt ${outpath}/esteff_phi${phi}_chr${chr}.txt
done
done
PRSCSSEL=${codepath}PRScs/PRScs_sel.R
Rscript ${PRSCSSEL} --plink_str ${PLINK} --valid_genotype ${val_geno} --valid_phenotype ${val_pheno} \
                    --esteff_path ${outpath} --index ${index} --outpath ${outpath}
fi

rm ${summary_file_prefix}.prscs

