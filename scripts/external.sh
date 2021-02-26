#!/bin/bash

while getopts ":s:e:p:o:" opt; do
  case $opt in
    s) summary_file="$OPTARG"
    ;;
    e) external_prefix="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    o) out_prefix="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument summary_file is %s  \033[0m\n" "$summary_file"
printf "\033[33mArgument external_prefix is %s  \033[0m\n" "$external_prefix"
printf "\033[33mArgument pop is %s  \033[0m\n" "$pop"
printf "\033[33mArgument out_prefix is %s  \033[0m\n" "$out_prefix"


EXTERNAL=/home/yasheng/comprsWeb/scripts/external.R
refLD=/home/yasheng/comprsWeb/1000GP/${pop}/LDmat_blk/

## External
Rscript ${EXTERNAL} --extsumm ${external_prefix} --esteff ${summary_file} --LDpath ${refLD} --outpath ${out_prefix}
