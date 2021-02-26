#!/bin/bash

DIR=/home/yasheng/comprsWeb/scripts/
SBLUP=${DIR}sblup.sh

# parameters
summary_file_prefix=/home/yasheng/comprsWeb/example_data/chr22/summary
ref_file=/home/yasheng/comprsWeb/example_data/chr22/geno
out_prefix=/home/yasheng/comprsWeb/example_data/output/SBLUP_esteff

sh ${SBLUP} -s ${summary_file_prefix} -H 0.1 -r ${ref_file} -t 1 -w 1000 -c 22 -o ${out_prefix}