#!/bin/bash

while getopts ":C:s:G:P:w:o:" opt; do
  case $opt in
    C) codepath="$OPTARG"
    ;;
    s) summary_file_prefix="$OPTARG"
    ;;
    G) val_geno="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    w) window_size="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument codepath is %s  \033[0m\n" "$codepath"
printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument valid_geno is %s  \033[0m\n" "$val_geno"
printf "\033[33mArgument valid_pheno is %s  \033[0m\n" "$val_pheno"
printf "\033[33mArgument p_len is %s  \033[0m\n" "$p_len"
printf "\033[33mArgument r2_val is %s  \033[0m\n" "$r2_val"
printf "\033[33mArgument dist_str is %s  \033[0m\n" "$dist_str"
if [ ${cov} -ne 0 ]
then 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"


# nps code
npsStr=${codepath}
plink=${codepath}plink
qctool=${codepath}qctool/qctool
cd ${npsStr}
dataID=ref

# window size
w1=`echo "${window_size}/4"|bc`
w2=`echo "${window_size}/2"|bc`
w3=`echo "${window_size}/4*3"|bc`

# val path
valchrpath=`echo ${val_geno%/*}`/
mkdir ${valchrpath}/tochr/

# tmp path
mkdir ${outpath}tmp

# Transform to nps format
## validation file (plink to dosage)
for chr in `seq 1 22`
do
mis=${valchrpath}tochr/chr${chr}
${plink} --silent --bfile ${val_geno} --chr ${chr} --make-bed --out ${mis}
bfile=${valchrpath}tochr/chr${chr}_imp
Rscript nps_imp.R --mis ${mis} --imp ${bfile}
${plink} --bfile ${bfile} --recode vcf --out ${bfile}
dosage=${valchrpath}tochr/chrom${chr}.${dataID}.dosage
${qctool} -g ${bfile}.vcf -filetype vcf -ofiletype dosage -og ${dosage} 
gzip -f ${dosage}
rm ${mis}*
done
# Standardize genotypes
./run_all_chroms.sh sge/nps_stdgt.job ${valchrpath}/tochr/ ${dataID}
cat ${valchrpath}/tochr/chrom*.${dataID}.snpinfo > ${valchrpath}/tochr/merge.${dataID}.snpinfo
## validation phenotype
Rscript nps_mk_ref.R --val ${val_geno} --valpheno ${val_pheno} --outpath ${valchrpath}tochr/
## summary statistics
Rscript nps_mk_summ.R --summ ${summary_file_prefix} --valpath ${valchrpath}/tochr/

# Initialize nps
Rscript npsR/nps_init.R --gwas ${summary_file_prefix}.nps\
               --train-dir ${valchrpath}/tochr/ --train-dataset ${dataID}\
               --window-size ${window_size} --out ${outpath}tmp --p 1

# Set up a special partition for GWAS-significant SNPs
echo step1
./run_all_chroms.sh sge/nps_gwassig.job ${outpath}tmp

# Set up the decorrelated "eigenlocus" space
echo step2
./run_all_chroms.sh sge/nps_decor_prune.job ${outpath}tmp 0
./run_all_chroms.sh sge/nps_decor_prune.job ${outpath}tmp ${w1}
./run_all_chroms.sh sge/nps_decor_prune.job ${outpath}tmp ${w2}
./run_all_chroms.sh sge/nps_decor_prune.job ${outpath}tmp ${w3}

# Partition the rest of genome
echo step3
Rscript npsR/nps_prep_part.R ${outpath}tmp 10 10
./run_all_chroms.sh sge/nps_part.job ${outpath}tmp 0
./run_all_chroms.sh sge/nps_part.job ${outpath}tmp ${w1}
./run_all_chroms.sh sge/nps_part.job ${outpath}tmp ${w2}
./run_all_chroms.sh sge/nps_part.job ${outpath}tmp ${w3}

# Estimate shrinkage weights for each partition
echo step4
Rscript npsR/nps_reweight.R ${outpath}tmp 1

# Summarize effect size
Rscript nps_est_eff.R --valpath ${valchrpath}/tochr/ --windowsize ${window_size} --tmppath ${outpath}tmp/ --outpath ${outpath}
