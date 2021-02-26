#!/bin/bash

CODEDIR=/home/yasheng/comprsWeb/scripts/
LDPRED2=${CODEDIR}LDpred2.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
summary_file_prefix=${DATADIR}chr22/summary
val_geno=${DATADIR}val/valid
val_pheno=${DATADIR}val/valid_pheno.txt

chr=22
outpath=${DATADIR}output/
sh ${LDPRED2} -s ${summary_file_prefix}.assoc.txt -G ${val_geno} -P ${val_pheno} -c 0 -C ${chr} -o ${outpath}
