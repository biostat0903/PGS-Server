#!/bin/bash

while getopts ":s:P:c:p:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    P) pop="$OPTARG"
    ;;
    c) chr="$OPTARG"
    ;;
    p) pi="$OPTARG"
    ;;
    o) out_prefix="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument ref_file is %s  \033[0m\n" "$ref_file"
printf "\033[33mArgument map_file is %s  \033[0m\n" "$map_file"
printf "\033[33mArgument chr is %s  \033[0m\n" "$chr"
printf "\033[33mArgument pi is %s  \033[0m\n" "$pi"
printf "\033[33mArgument out_prefix is %s  \033[0m\n" "$out_prefix"

CURRENT_DIR=$(cd `dirname $0`; pwd)
PACKAGE_DIR=/home/yasheng/comprsWeb/software/
GCTB=${PACKAGE_DIR}/gctb_2.0_Linux/gctb
refLD=/home/yasheng/comprsWeb/1000GP/${pop}/LDmat_shrinkage/chr${chr}

## SbayesR
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summary_file_prefix}.assoc.txt > ${summary_file_prefix}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summary_file_prefix}.ma

${GCTB} --sbayes R --ldm ${refLD}.ldm.shrunk --pi ${pi} \
        --gamma 0.0,0.01,0.1,1 --gwas-summary ${summary_file_prefix}.ma \
        --exclude-mhc --chain-length 10000 --burn-in 2000 --out-freq 5000 --out ${out_prefix}

## remove file
rm ${summary_file_prefix}.ma
rm ${out_prefix}.parRes
rm ${out_prefix}.mcmcsamples.SnpEffects
rm ${out_prefix}.mcmcsamples.Par
rm ${out_prefix}.mcmcsamples.CovEffects
rm ${out_prefix}.covRes
