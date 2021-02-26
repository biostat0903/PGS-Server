#!/bin/bash

CODEDIR=/home/yasheng/comprsWeb/scripts/
CT=${CODEDIR}CT.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
summ=${DATADIR}all/summary.assoc.txt
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
pl=2
rv=0.1,0.2
dv=50
outpath=${DATADIR}output/
sh ${CT} -s ${summ} -G ${valg} -P ${valp} -p ${pl} -r ${rv} -d ${dv} -c 0 -t 1 -o ${outpath}

