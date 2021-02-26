#!/bin/bash

DIR=/home/yasheng/comprsWeb/scripts/
DBSLMM=${DIR}DBSLMM1.sh

# parameters
DATADIR=/home/yasheng/comprsWeb/example_data/
summ=${DATADIR}all/summary
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
herit=0.1
type=d
cov=0
index=r2
thread=1
outpath=/home/yasheng/comprsWeb/example_data/output/
sh ${DBSLMM} -s ${summ} -H ${herit} -m 1062150 -n 300 -p EUR\
 -G ${valg} -P ${valp} -T ${type} -c ${cov} -i ${index} -t ${thread} -o ${outpath}
