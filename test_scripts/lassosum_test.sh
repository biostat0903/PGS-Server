#!/bin/bash

CODEDIR=/home/yasheng/comprsWeb/scripts/
LASSOSUM=${CODEDIR}lassosum.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
summ=${DATADIR}all/summary
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
pop=EUR

sh ${LASSOSUM} -s ${summ} -G ${valg} -P ${valp} -p ${pop} -c 0 -o ${outpath}

