# COM-PGS
The software compares the existing methods of PGS construction, including CT (`bigsnpr`), DBSLMM, lassosum, LDpred2-auto, LDpred2-inf, LDpred2-sp, LDpred2-nosp, NPS, PRSCS, SbayesR, SBLUP and SCT. Specially, for DBSLMM and PRSCS, we give two versions: automatic and tuning version. The code for each method is in the `script` file folder. 

## Preparing for code
All methods are coded by `R`, `plink` and shell script. The user should run the code in Linux or virtual environment for Linux. Then the user should install the R packages, including `bigsnpr`, `bigstatsr`, `bigreadr`, `plyr`, `tidyverse`, `optparse`, `lassosum` and `doParallel`. The user also need to download the [plink] (https://www.cog-genomics.org/plink/). 
## CT (`bingsnpr`)
We use the recommendation setting of `bigsnpr`, including 1,400 parameter combinations. The combination has three parameters, including p value, window size and r2. follwing , we recomment 50 different settings of p value, 4 different settings of window size and 7 different settings of r2. 
The script `CT.sh` is to call `CT.R` function. The shell script is as following:
````shell
# code path
CODEDIR=/home/yasheng/comprsWeb/scripts/
DATADIR=/home/yasheng/comprsWeb/example_data/
CT=${CODEDIR}CT.sh
# data path
summ=${DATADIR}all/summary.assoc.txt
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# paramaters
pl=2
rv=0.1,0.2
dv=50
# CT method
sh ${CT} -C ${CODEDIR}01_CT/ -s ${summ} -G ${valg} -P ${valp} -p ${pl} -r ${rv} -d ${dv} -${outpath}
````

## DBSLMM (`R`+`C++`+`shell`)
The detail for DBSLMM is: https://github.com/biostat0903/DBSLMM. 

## lassosum (`lassosum` R package)
We use the default setting for lassosum using `R`. Lassosum contains two hyper-parameters: the penalty parameter (λ) in the lasso regression and the shrinkage parameter (s) used for computing the SNP correlation matrix in the reference panel. Following lassosum paper, we considered four choices of s (0.2, 0.5, 0.9 and 1) and 20 choices of λ that are evenly spaced on the logarithmic scale between log(0.01) and log(0.1). The script `lassosum.sh` is to call `lassosum.R` function. The shell script is as following:
````shell
# code path
CODEDIR=/home/yasheng/comprsWeb/scripts/
LASSOSUM=${CODEDIR}lassosum.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
# data path
summ=${DATADIR}all/summary
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# population
pop=EUR
# lassosum method
sh ${LASSOSUM} -C ${CODEDIR}/03_lassosum -s ${summ} -G ${valg} -P ${valp} -p ${pop} -o ${outpath}
````

## LDpred2 (`bigsnpr`)
Following LDpred2 paper, we examined four different models implemented in LDpred2 described as follows. (1) LDpred2-inf is the infinitesimal model that is fitted based on an analytic solution. (2) LDpred2-sp is a sparse Bayesian variable selection regression model that selects a small proportion of SNPs to construct PGS. LDpred2-sp contains two hyper-parameters that include the proportion of causal variants p and the SNP heritability h2. LDpred2-sp explores different combinations of the two hyper-parameters on a pre-selected set of grid values and determines the optimal hyper-parameter combination through cross-validation. (3) LDpred2-nosp fits the same model as LDpred2-sp but sets the proportion of causal variants p to be exactly one (and thus becomes non-sparse). (4) LDpred2-auto fits the same model as LDpred2-nosp but automatically estimates p and h2 from the training data. The script `LDpred2.sh` is to call `LDpred2.R` function. The shell script is as following:
````shell
# code path
CODEDIR=/home/yasheng/comprsWeb/scripts/
LDPRED2=${CODEDIR}LDpred2.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
# data path
summary_file_prefix=${DATADIR}chr22/summary
val_geno=${DATADIR}val/valid
val_pheno=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# parameters
chr=22
# LDpred2 method
sh ${LDPRED2} -C ${CODEDIR}/04_LDpred2 -s ${summary_file_prefix}.assoc.txt -G ${val_geno} -P ${val_pheno} -C ${chr} -o ${outpath}
````

## NPS (`C++`+`R`)
The original version of `NPS` can not use to analyze multiple traits at one time and can not analyze the genotype with NA. We update the `NPS` pacakge. The script `LDpred2.sh` is to call `LDpred2.R` function. The shell script is as following:
````shell
# code path
CODEDIR=/home/yasheng/comprsWeb/scripts/
NPS=${CODEDIR}nps.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
# data path
val_geno=/home/yasheng/comprsWeb/example_data/val/valid
val_pheno=/home/yasheng/comprsWeb/example_data/val/valid_pheno.txt
summary_file_prefix=/home/yasheng/comprsWeb/example_data/all/summary
outpath=/home/yasheng/comprsWeb/example_data/output/
# parameter
window_size=60
# NPS
sh ${nps} -C ${CODEDIR}/05_NPS -s ${summary_file_prefix} -G ${val_geno} -P ${val_pheno} -w ${window_size} -o ${outpath}
````

## PRSCS
Following PRSCS paper, we set the hyper-parameter a in PRSCS to the default value of 1, set the hyper-parameter b to the default value of 0.5, and inferred the global scaling hyper-parameter ϕ among a set of four choices {10^(-6),10^(-4),0.01,1}. We also examined the automatic version of PRSCS, referring to PRSCS-auto. 
````shell

DIR=/home/yasheng/comprsWeb/scripts/
PRSCS=${DIR}PRSCS.sh

# parameters
summary_file_prefix=/home/yasheng/comprsWeb/example_data/chr22/summary
out_prefix=/home/yasheng/comprsWeb/example_data/output/PRSCS_esteff
chr=22
pop=EUR

sh ${PRSCS} -s ${summary_file_prefix} -c ${chr} -p ${pop} -o ${out_prefix}


````

## SbayesR
````bash
# code path
CODEDIR=/home/yasheng/comprsWeb/scripts/07_SbayesR/
SBAYESR=${DIR}SbayesR.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
# parameters
summary_file_prefix=${DATADIR}chr22/summary
chr=22
pi=0.95,0.02,0.02,0.01
out_prefix=${DATADIR}output/SbayesR_esteff
pop=EUR

# SbayesR method
sh ${SBAYESR} -C ${CODEDIR} -s ${summary_file_prefix} -P ${pop} -c ${chr} -p ${pi} -o ${out_prefix}
````
## SBLUP

## SCT
