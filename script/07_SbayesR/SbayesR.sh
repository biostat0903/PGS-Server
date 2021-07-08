#!/bin/bash

while getopts ":C:s:L:c:p:o:" opt; do
  case $opt in
    C) codepath="$OPTARG"
    ;;
    s) summary_file_prefix="$OPTARG"
    ;;
    L) LDpath="$OPTARG"
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

printf "\033[33mArgument codepath is %s  \033[0m\n" "$codepath"
printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument LDpath is %s  \033[0m\n" "$LDpath"
printf "\033[33mArgument chr is %s  \033[0m\n" "$chr"
printf "\033[33mArgument pi is %s  \033[0m\n" "$pi"
printf "\033[33mArgument out_prefix is %s  \033[0m\n" "$out_prefix"

GCTB=${codepath}/gctb_2.0_Linux/gctb
refLD=${LDpath}/chr${chr}

## SbayesR
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summary_file_prefix}.assoc.txt > ${summary_file_prefix}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summary_file_prefix}.ma

outputPath=`echo ${out_prefix%/*}`
${GCTB} --sbayes R --ldm ${refLD}.ldm.shrunk --pi ${pi} \
        --gamma 0.0,0.01,0.1,1 --gwas-summary ${summary_file_prefix}.ma \
        --exclude-mhc --chain-length 10000 --burn-in 2000 --out-freq 5000 --out ${outputPath}/summary
		
mv  ${outputPath}/summary.snpRes ${outputPath}/sbayesr_esteff.txt
## remove file
rm ${summary_file_prefix}.ma
rm ${outputPath}/summary.parRes
rm ${outputPath}/summary.mcmcsamples.SnpEffects
rm ${outputPath}/summary.mcmcsamples.Par
rm ${outputPath}/summary.mcmcsamples.CovEffects
rm ${outputPath}/summary.covRes
