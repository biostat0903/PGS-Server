#!/bin/bash

while getopts "s:H:m:n:p:G:P:T:c:i:t:o:" opt; do
  case $opt in
    s) summary_file_prefix="$OPTARG"
    ;;
    H) herit="$OPTARG"
    ;;
    m) nsnp="$OPTARG"
    ;;
    n) n="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    G) val_geno_prefix="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    T) type="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    i) index="$OPTARG"
    ;;
    t) thread="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument herit is %s  \033[0m\n" "$herit"
printf "\033[33mArgument val_geno_prefix is %s  \033[0m\n" "$val_geno_prefix"
printf "\033[33mArgument valid_pheno is %s  \033[0m\n" "$val_pheno"
printf "\033[33mArgument type is %s  \033[0m\n" "$type"
if [ ${cov} -ne 0 ]; then 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument index is %s  \033[0m\n" "$index"
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"

SPLIT=/home/yasheng/comprsWeb/software/DBSLMM/SPLITCHR.R
DBSLMM=/home/yasheng/comprsWeb/software/DBSLMM/DBSLMM.R
TUNE=/home/yasheng/comprsWeb/software/DBSLMM/TUNE.R
plink=/home/yasheng/comprsWeb/software/plink
dbslmm=/home/yasheng/comprsWeb/software/DBSLMM/dbslmm
block_prefix=/home/yasheng/comprsWeb/block_data/${pop}/chr
val_geno=${val_geno_prefix}

# split to chromosome
Rscript ${SPLIT} --summary ${summary_file_prefix}.assoc.txt

# DBSLMM: tuning version
if [ "$type" = "t" ]
then

for chr in `seq 1 22`
do

BLOCK=${block_prefix}${chr}
summchr=${summary_file_prefix}_chr${chr}
Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${plink}\
                  --dbslmm ${dbslmm} --ref ${val_geno} --n ${n} --type ${type} --nsnp ${nsnp} --block ${BLOCK}.bed\
                  --h2 ${herit} --h2f 0.8,1,1.2 --thread ${thread}
summchr_prefix=`echo ${summchr##*/}`
for h2f in 0.8 1 1.2
do
if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps" ];then
rm ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps
fi
${plink}  --silent --bfile ${val_geno} --score ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.txt 1 2 4 sum\
          --out ${outpath}${summchr_prefix}_h2f${h2f}
rm ${outpath}${summchr_prefix}_h2f${h2f}.log
if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.nopred" ];then
rm ${outpath}${summchr_prefix}_h2f${h2f}.nopred
fi
done

done

summchr_prefix2=`echo ${summchr_prefix%_*}`
if [ ${cov} -eq 0 ]
then 
Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
       --h2Range 0.8,1,1.2 --index ${index}
else 
Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
       --h2Range 0.8,1,1.2 --index ${index} --cov ${cov}
fi

hbest=`cat ${outpath}${summchr_prefix2}_hbest.${index}`
for chr in `seq 1 22`
do
mv ${outpath}${summchr_prefix2}_chr${chr}_h2f${hbest}.dbslmm.txt ${outpath}${summchr_prefix2}_chr${chr}_best.dbslmm.txt
rm ${outpath}${summchr_prefix2}_chr${chr}_h2f*
done

fi

# DBSLMM default version
if [ "$type" = "d" ]
then
for chr in `seq 1 22` 
do
BLOCK=${block_prefix}${chr}
summchr=${summary_file_prefix}_chr${chr}

Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${plink}\
                  --dbslmm ${dbslmm} --ref ${val_geno} --n ${n} --nsnp ${nsnp} --block ${BLOCK}.bed\
                  --h2 ${herit} --thread ${thread}
summchr_prefix=`echo ${summchr##*/}`
rm ${outpath}${summchr_prefix}.dbslmm.badsnps

done
fi
