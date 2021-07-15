# COM-PGS
The software compares the existing methods of PGS construction, including CT (`bigsnpr`), DBSLMM, lassosum, LDpred2-auto, LDpred2-inf, LDpred2-sp, LDpred2-nosp, NPS, PRSCS, SbayesR, SBLUP and SCT. Specially, for DBSLMM and PRSCS, we give two versions: automatic and tuning version. The code for each method is in the `script` file folder. 

## Preparing for code
All methods are coded by `R`, `plink` and shell script. The user should run the code in Linux or virtual environment for Linux. Then the user should install the R packages, including `bigsnpr`, `bigstatsr`, `bigreadr`, `plyr`, `tidyverse`, `optparse`, `lassosum` and `doParallel`. <br>
NPS should use [`qctool`](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/download.html) and [`plink`](https://www.cog-genomics.org/plink/). DBSLMM should use `plink`. SBayesR and SBLUP should use [`GCTB`](https://cnsgenomics.com/software/gctb/#Download) and [`GCTA`](https://cnsgenomics.com/software/gcta/#Download). 
## CT (bingsnpr `R` package)
We use the recommendation setting of `bigsnpr`, including 1,400 parameter combinations. The combination has three parameters, including p value, window size and r2. follwing , we recomment 50 different settings of p value, 4 different settings of window size and 7 different settings of r2. 
The script `CT.sh` is to call `CT.R` function. The shell script is as following:
````shell
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}script/
DATADIR=${DIR}/example_data/
CT=${CODEDIR}01_CT/CT.sh
# data path
summ=${DATADIR}all/summary
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# paramaters
pl=2
rv=0.1,0.2
dv=50
# CT method
sh ${CT} -C ${CODEDIR}01_CT/ -s ${summ} -G ${valg} -P ${valp} -p ${pl} -r ${rv} -d ${dv} -o ${outpath}
````

## DBSLMM (`R`+`C++`+`shell`)
The detail for DBSLMM is: https://github.com/biostat0903/DBSLMM. 

## lassosum (lassosum `R` package)
We use the default setting for lassosum using `R`. Lassosum contains two hyper-parameters: the penalty parameter (λ) in the lasso regression and the shrinkage parameter (s) used for computing the SNP correlation matrix in the reference panel. Following lassosum paper, we considered four choices of s (0.2, 0.5, 0.9 and 1) and 20 choices of λ that are evenly spaced on the logarithmic scale between log(0.01) and log(0.1). The script `lassosum.sh` is to call `lassosum.R` function. The shell script is as following:
````shell
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}script/
DATADIR=${DIR}/example_data/
LASSOSUM=${CODEDIR}03_lassosum/lassosum.sh
# data path
summ=${DATADIR}all/summary
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# population
pop=EUR
# lassosum method
sh ${LASSOSUM} -C ${CODEDIR}03_lassosum/ -s ${summ} -G ${valg} -P ${valp} -p ${pop} -o ${outpath}
````

## LDpred2 (bingsnpr `R` package)
Following LDpred2 paper, we examined four different models implemented in LDpred2 described as follows. (1) LDpred2-inf is the infinitesimal model that is fitted based on an analytic solution. (2) LDpred2-sp is a sparse Bayesian variable selection regression model that selects a small proportion of SNPs to construct PGS. LDpred2-sp contains two hyper-parameters that include the proportion of causal variants p and the SNP heritability h2. LDpred2-sp explores different combinations of the two hyper-parameters on a pre-selected set of grid values and determines the optimal hyper-parameter combination through cross-validation. (3) LDpred2-nosp fits the same model as LDpred2-sp but sets the proportion of causal variants p to be exactly one (and thus becomes non-sparse). (4) LDpred2-auto fits the same model as LDpred2-nosp but automatically estimates p and h2 from the training data. The script `LDpred2.sh` is to call `LDpred2.R` function. The shell script is as following:
````shell
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}script/
DATADIR=${DIR}/example_data/
LDPRED2=${CODEDIR}/04_LDpred2/LDpred2.sh
# data path
summary_file_prefix=${DATADIR}chr22/summary
val_geno=${DATADIR}val/valid
val_pheno=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# parameters
chr=22
# LDpred2 method
sh ${LDPRED2} -C ${CODEDIR}04_LDpred2/ -s ${summary_file_prefix} -G ${val_geno} -P ${val_pheno} -R ${chr} -o ${outpath}
````

## NPS (`C++`+`R`)
The original version of `NPS` can not use to analyze multiple traits at one time and can not analyze the genotype with NA. We update the `NPS` pacakge. The script `LDpred2.sh` is to call `LDpred2.R` function. The shell script is as following:
````shell
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}script/
DATADIR=${DIR}/example_data/
NPS=${CODEDIR}05_NPS/NPS.sh
SOFTDIR=${DIR}/software/
# data path
val_geno=${DATADIR}/val/valid
val_pheno=${DATADIR}/val/valid_pheno.txt
summary_file_prefix=${DATADIR}/all/summary
outpath=${DATADIR}/output/
# parameter
window_size=60
# NPS
sh ${NPS} -C ${SOFTDIR} -s ${summary_file_prefix} -G ${val_geno} -P ${val_pheno} -w ${window_size} -o ${outpath}
````

## PRSCS (`python`)
Following PRSCS paper, we set the hyper-parameter a in PRSCS to the default value of 1, set the hyper-parameter b to the default value of 0.5, and inferred the global scaling hyper-parameter ϕ among a set of four choices {10^(-6),10^(-4),0.01,1}. We also examined the automatic version of PRSCS, referring to PRSCS-auto. 
````shell
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}software/
DATADIR=${DIR}/example_data/
PRSCS=${DIR}/script/06_PRSCS/PRSCS.sh

# parameters
summary_file_prefix=${DATADIR}all/summary
mkdir ${DATADIR}/output/PRSCS
outpath=${DATADIR}/output/PRSCS
index=r2
LDpath=${DATADIR}all/ldblk_1kg_eur
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt

# PRSCS-auto
sh ${PRSCS} -C ${CODEDIR} -s ${summary_file_prefix} -L ${LDpath} -T auto \
            -G ${valg} -o ${outpath}
# PRSCS-tuning
sh ${PRSCS} -C ${CODEDIR} -s ${summary_file_prefix} -L ${LDpath} -T tuning \
            -G ${valg} -o ${outpath} -P ${valp} -i ${index}
````
Note: The user should download the LD matrix on [PRSCS](https://github.com/getian107/PRScs).

## SbayesR (GCTB `C++`)
Following SbayesR paper, we set the weights of the four normal components (“--pi”) to the default values of {0.95,0.02,0.02,0.01} and set the four scaling variance parameters (“--gamma”) to the default values of {0,0.01,0.1,1}. We constructed the SNP LD matrix using the “--make-shrunk-ldm” option, again with the default settings (effective population size = 11,400; genetic map sample size = 183; shrinkage hard threshold = 10-5). We set the MCMC chain length to be 10,000 with an additional 2,000 burn-in iterations. 
````bash
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}script/
DATADIR=${DIR}/example_data/
SBAYESR=${CODEDIR}07_SbayesR/SbayesR.sh
SOFTDIR=${DIR}/software/ # please download the GCTB software and store it in the directory
LDDIR=${DATADIR}chr22
# parameters
summary_file_prefix=${DATADIR}chr22/summary
chr=22
pi=0.95,0.02,0.02,0.01
out_prefix=${DATADIR}output/SbayesR_esteff
pop=EUR
# SbayesR method
sh ${SBAYESR} -C ${SOFTDIR} -s ${summary_file_prefix} -L ${LDDIR} -c ${chr} -p ${pi} -o ${out_prefix}
````
Note: The user should make the LD matrix following the manual of [SbayesR](https://cnsgenomics.com/software/gctb/#Tutorial).

## SBLUP (GCTA `C++`)
We used the GCTA to fit SBLUP and used h2 as the SNP heritability input. SBLUP also requires users to specify a LD window size that is used to construct the SNP correlation matrix in the reference panel. 
````bash
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}/script/
DATADIR=${DIR}/example_data/
SBLUP=${CODEDIR}/08_SBLUP/SBLUP.sh
SOFTDIR=${DIR}/software/ # please download the GCTA software and store it in the directory
# parameters
summary_file_prefix=${DATADIR}chr22/summary
herit=0.1
window=1000
chr=22
ref_file=${DATADIR}chr22/geno
out_prefix=${DATADIR}output/SBLUP_esteff
# SBLUP method
sh ${SBLUP} -C ${SOFTDIR} -s ${summary_file_prefix} -H ${herit} -r ${ref_file} -t 1 -w ${window} -c ${chr} -o ${out_prefix}
````

## SCT (bingsnpr `R` package)
The input of SCT is the same as that of CT.
````shell
# code path
DIR=/public/home/biostat03/project/compProject/PGS-Server-main/
CODEDIR=${DIR}script/
DATADIR=${DIR}/example_data/
SCT=${CODEDIR}09_SCT/SCT.sh
# data path
summ=${DATADIR}all/summary
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
outpath=${DATADIR}output/
# paramaters
pl=2
rv=0.1,0.2
dv=50
# SCT method
sh ${SCT} -C ${CODEDIR}09_SCT/ -s ${summ} -G ${valg} -P ${valp} -p ${pl} -r ${rv} -d ${dv} -o ${outpath}
````
## External validation (`R`)
The detail for external validation is: https://github.com/biostat0903/DBSLMM. 

