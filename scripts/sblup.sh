#!/bin/bash

while getopts ":s:H:r:t:w:c:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    H) heritability="$OPTARG"
    ;;
    r) ref_file="$OPTARG"
    ;;
    t) thread="$OPTARG"
    ;;
    w) window_size="$OPTARG"
    ;;
    c) chr="$OPTARG"
    ;;
    o) out_prefix="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument heritability is %s  \033[0m\n" "$heritability"
printf "\033[33mArgument ref_file is %s  \033[0m\n" "$ref_file"
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
printf "\033[33mArgument window_size is %s  \033[0m\n" "$window_size"
printf "\033[33mArgument out_prefix is %s  \033[0m\n" "$out_prefix"

CURRENT_DIR=$(cd `dirname $0`; pwd)
PACKAGE_DIR=/home/yasheng/comprsWeb/software/

## snp number
m=`cat ${ref_file}.bim | wc -l`
echo "m: $m"
cojo=$(echo "${m}*(1/${heritability}-1)" | bc -l)
echo "cojo: $cojo"

## sblup estimation
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summary_file_prefix}.assoc.txt > ${summary_file_prefix}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summary_file_prefix}.ma

${PACKAGE_DIR}gcta_1.93.2beta/gcta64 --bfile ${ref_file} --chr ${chr} --cojo-file ${summary_file_prefix}.ma --cojo-sblup ${cojo} --cojo-wind 1000 --thread-num 1 --out ${out_prefix}

## remove file
rm ${out_prefix}.*badsnps
rm ${summary_file_prefix}.ma
rm ${out_prefix}.log

exit 0
