#!/bin/bash
#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=dos
#SBATCH --mem=10G
#SBATCH --cpus-per-task=6

#SBATCH --array=1-21
#SBATCH --output=/net/mulan/disk2/yasheng/comparisonProject/00_cluster_file/09_todosage_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/comparisonProject/00_cluster_file/09_todosage_%a.err



# dataID=ref
# qctool=/net/mulan/home/yasheng/comparisonProject/program/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool
# for chr in `seq 1 21`
# do 
# bfile=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/geno/chr${chr}_imp
# dosage=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage/chrom${chr}.${dataID}.dosage
# plink-1.9 --bfile ${bfile} --recode vcf --out ${bfile}
# ${qctool} -g ${bfile}.vcf -filetype vcf -ofiletype dosage -og ${dosage}
# gzip ${dosage}
# rm ${bfile}.vcf
# done

# ## Standardize genotypes
# refpath=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage
# cd /net/mulan/disk2/yasheng/comparisonProject/nps-master/
# ./run_all_chroms.sh sge/nps_stdgt.job ${refpath} ${dataID}
# cd /net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage/
# cat chrom*.${dataID}.snpinfo > merge.${dataID}.snpinfo


dataID=val
qctool=/net/mulan/home/yasheng/comparisonProject/program/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool
imp=/net/mulan/disk2/yasheng/comparisonProject/code/02_method/09_nps_imp.R

SEND_THREAD_NUM=21
tmp_fifofile="/tmp/$$.fifo"
mkfifo "$tmp_fifofile"
exec 6<>"$tmp_fifofile"

for ((i=0;i<$SEND_THREAD_NUM;i++));do
echo
done >&6

for chr in `seq 1 21`
do 
read -u6
{

## snp id
refsnp=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage/chrom${chr}.ref.snpinfo
refid=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage/chrom${chr}.id.txt
awk '{print $2}' ${refsnp} > ${refid}
sed -i '1d' ${refid}

## select snp
mis1=/net/mulan/disk2/yasheng/predictionProject/plink_file/hm3/chr${chr}
mis2=/net/mulan/disk2/yasheng/predictionProject/plink_file/hm3/xchr${chr}
plink-1.9 --bfile ${mis1} --extract ${refid} --make-bed --out ${mis2}
bfile=/net/mulan/disk2/yasheng/predictionProject/plink_file/hm3/chr${chr}_imp
rm ${mis2}.bk
Rscript ${imp} --mis ${mis2} --imp ${bfile}
dosage=/net/mulan/disk2/yasheng/predictionProject/plink_file/dosage/chrom${chr}.${dataID}.dosage
plink-1.9 --bfile ${bfile} --recode vcf --out ${bfile}
${qctool} -g ${bfile}.vcf -filetype vcf -ofiletype dosage -og ${dosage} 
gzip -f ${dosage}

## remove file
# rm ${bfile}.vcf
# rm ${bfile}.bed
# rm ${bfile}.bim
# rm ${bfile}.fam
rm ${mis}.bk
rm ${mis}.rds

} &
pid=$!
echo $pid
done

wait

exec 6>&-

qctool=/net/mulan/home/yasheng/comparisonProject/program/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool

chr=22
refsnp=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage/chrom${chr}.ref.snpinfo
refid=/net/mulan/disk2/yasheng/comparisonProject/03_subsample/hm3/dosage/chrom${chr}.id.txt
awk '{print $2}' ${refsnp} > ${refid}
sed -i '1d' ${refid}
dosage1=/net/mulan/disk2/yasheng/predictionProject/plink_file/dosage/chrom${chr}.val.dosage
dosage2=/net/mulan/disk2/yasheng/predictionProject/plink_file/dosage/chrom${chr}.valid.dosage
gunzip ${dosage1}.gz
${qctool} -g ${dosage1} -og ${dosage2} -incl-variants ${refid}
gzip ${dosage2}


