# COM-PGS
The software compares the existing methods of PGS construction, including CT (`bigsnpr`), DBSLMM, lassosum, LDpred2-auto, LDpred2-inf, LDpred2-sp, LDpred2-nosp, NPS, PRSCS, SbayesR, SBLUP and SCT. Specially, for DBSLMM and PRSCS, we give two versions: automatic and tuning version. The code for each method is in the `script` file folder. 
## Preparing for code
All methods are coded by `R`, `plink` and shell script. The user should run the code in Linux or virtual environment for Linux. Then the user should install the R packages, including `bigsnpr`, `bigstatsr`, `bigreadr`, `plyr`, `tidyverse`, `optparse`, `lassosum` and `doParallel`. The user also need to download the [plink] (https://www.cog-genomics.org/plink/). 
## CT (`bingsnpr`)
We use the recommendation setting of `bigsnpr`, including 1,400 parameter combinations. The combination has three parameters, including p value, window size and r2. follwing , we recomment 50 different settings of p value, 4 different settings of window size and 7 different settings of r2. 
The script `CT.sh` is to call `CT.R` function. The shell script is as following:
````shell
CODEDIR=/home/yasheng/comprsWeb/scripts/
CT=${CODEDIR}CT.sh
DATADIR=/home/yasheng/comprsWeb/example_data/
summ=${DATADIR}all/summary.assoc.txt
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
pl=2
rv=0.1,0.2
dv=50
outpath=${DATADIR}output/
sh ${CT} -C ${CODEDIR} -s ${summ} -G ${valg} -P ${valp} -p ${pl} -r ${rv} -d ${dv} -${outpath}
````
