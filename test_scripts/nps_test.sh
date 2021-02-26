#!/bin/bash
nps=/home/yasheng/comprsWeb/scripts/nps.sh

val_geno=/home/yasheng/comprsWeb/example_data/val/valid
val_pheno=/home/yasheng/comprsWeb/example_data/val/valid_pheno.txt
summary_file_prefix=/home/yasheng/comprsWeb/example_data/all/summary
window_size=60
outpath=/home/yasheng/comprsWeb/example_data/output/

sh ${nps} -s ${summary_file_prefix} -G ${val_geno} -P ${val_pheno} -w ${window_size} -o ${outpath}