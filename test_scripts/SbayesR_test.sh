#!/bin/bash

DIR=/home/yasheng/comprsWeb/scripts/
SBAYESR=${DIR}SbayesR.sh

# parameters
summary_file_prefix=/home/yasheng/comprsWeb/example_data/chr22/summary
chr=22
pi=0.95,0.02,0.02,0.01
out_prefix=/home/yasheng/comprsWeb/example_data/output/SbayesR_esteff
pop=AFR
sh ${SBAYESR} -s ${summary_file_prefix} -P ${pop} -c ${chr} -p ${pi} -o ${out_prefix}
