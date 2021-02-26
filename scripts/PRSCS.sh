#!/bin/bash

while getopts ":s:c:p:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    c) chr="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    o) out_prefix="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument chr is %s  \033[0m\n" "$chr"
printf "\033[33mArgument population is %s  \033[0m\n" "$pop"
printf "\033[33mArgument out_prefix is %s  \033[0m\n" "$out_prefix"

CURRENT_DIR=$(cd `dirname $0`; pwd)
PACKAGE_DIR=/home/yasheng/comprsWeb/software/PRScs/
PRSCS=${PACKAGE_DIR}/PRScs.py
popl=`echo $pop | tr 'A-Z' 'a-z'`
REF=/home/yasheng/comprsWeb/1000GP/${pop}/ldblk_1kg_${popl}

awk '{print $2,$6,$7,$9,$11}' ${summary_file_prefix}.assoc.txt > ${summary_file_prefix}.prscs
sed -i '1i\SNP A1 A2 BETA P' ${summary_file_prefix}.prscs
nobs=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
bimStr=/home/yasheng/comprsWeb/1000GP/${pop}/geno/chr${chr}
python3 ${PRSCS} --ref_dir=${REF} --bim_prefix=${bimStr} --sst_file=${summary_file_prefix}.prscs --n_gwas=${n} --chrom=${chr} --out_dir=${out_prefix}
mv ${out_prefix}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt ${out_prefix}PRSCS_esteff_chr${chr}.txt 
rm ${summary_file_prefix}.prscs
