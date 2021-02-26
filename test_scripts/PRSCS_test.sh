#!/bin/bash

DIR=/home/yasheng/comprsWeb/scripts/
PRSCS=${DIR}PRSCS.sh

# parameters
summary_file_prefix=/home/yasheng/comprsWeb/example_data/chr22/summary
out_prefix=/home/yasheng/comprsWeb/example_data/output/PRSCS_esteff
chr=22
pop=EUR

sh ${PRSCS} -s ${summary_file_prefix} -c ${chr} -p ${pop} -o ${out_prefix}
