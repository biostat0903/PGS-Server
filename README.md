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
sh ${LASSOSUM} -C ${CODEDIR}/03_lassosum -s ${summ} -G ${valg} -P ${valp} -p ${pop} -c 0 -o ${outpath}
````

## LDpred2 (`bigsnpr`)
