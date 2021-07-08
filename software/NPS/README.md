# Non-Parametric Shrinkage (NPS)
NPS implements a non-parametric polygenic risk prediction algorithm described in [Chun et al. 2020 (preprint)](https://doi.org/10.1101/370064). NPS transforms genetic data into an orthogonal domain called "eigenlocus space". Then, it re-weights GWAS effect sizes by partitioning genetic variations into trenches and measuring the predictive power of each trench in an independent training cohort. To run NPS, two sets of data are required: GWAS summary statistics and small individual-level training cohort with both genotype and phenotype data. 

We recommend running NPS in computer clusters. The following platforms are supported: 
* SGE/UGER
* LSF 
* Slurm
* MacOS desktop (not recommended for genome-wide datasets)
* Linux desktop (not recommended for genome-wide datasets)

For inquiries on software, please contact: 
* Sung Chun (SungGook.Chun@childrens.harvard.edu)
* Nathan Stitziel (nstitziel@wustl.edu) 
* Shamil Sunyaev (ssunyaev@rics.bwh.harvard.edu). 

## How to Install
1. Download and unpack NPS package ([version 1.1.1](https://github.com/sgchun/nps/archive/1.1.1.tar.gz)). Part of NPS codes are optimized in C++ and have to be compiled using GNU C++ compiler (GCC-4.4 or higher). This will create two executable binaries, **stdgt** and **grs**, in the top-level directory. Note: On Mac OS X, [Command Line Tools for Xcode](https://developer.apple.com/download/more/) has to be installed first in order to use gcc and make.
   ```bash
   tar -zxvf nps-1.1.1.tar.gz
   cd nps-1.1.1/
   make
   ```

2. The core NPS module was implemented in R (version 3.3 or higher). Although NPS can run on a standard version of R, we strongly recommend using R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [Microsoft R open](https://mran.microsoft.com/open). These libraries can substantially speed up NPS operations.  

3. NPS relies on R libraries, [pROC](https://cran.r-project.org/web/packages/pROC/index.html) and [DescTools](https://cran.r-project.org/web/packages/DescTools/index.html), to report the accuracy of polygenic scores in AUC and Nagelkerke's R^2. These modules are optional; if they are not available, AUC and Nagelkerke's R^2 calculation will be skipped. To install these packages, run the following on command line: 
   ```bash
   Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
   Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
   ```

   To install the R extensions in the home directory (e.g. ~/R) rather than in the default system path, use the following commands instead:
   ```bash
   Rscript -e 'install.packages("pROC", "~/R", repos="http://cran.r-project.org")' 
   Rscript -e 'install.packages("DescTools", "~/R", repos="http://cran.r-project.org")' 
   ```
   Then, add "~/R" to the local R library path in your login shell's start-up file. For example, in case of bash, add the following to .bash_profile or .bashrc: 
   ```bash
   export R_LIBS="~/R:$R_LIBS"
   ```

## NPS Test datsets
We provide two sets of simulated test cases. They are provided separately from the software distribution and can be downloaded from Sunyaev Lab FTP server (**ftp://genetics.bwh.harvard.edu/download/schun/**). Test set #1 is relatively small and can be easily run on desktop whereas test set #2 is a more realistic dataset but requires serious computational resource to run NPS.

Test set | Total # of simulated SNPs | # of simulated causal SNPs | NPS disk space requirement | NPS run time
--- | --- | --- | --- | --- 
#1 | 100,449 | 522 (0.5%) | 5 GB | 30 mins* 
#2 | 5,012,500 | 5,008 (0.1%) | 800 GB | 3-6 hours**

[*] On a desktop computer, without parallelization  
[**] On computer clusters, parallelizing on up to 88 CPUs, with linear algebra acceleration.  

Download the test dataset and unpacked it as below. See [File Formats](https://github.com/sgchun/nps/blob/master/FileFormats.md) for the description on the input file formats.
```bash
cd nps-1.1.1/testdata/

tar -xvf NPS.Test1.tar  
# This will create the following test data files in nps-1.1.1/testdata/Test1
# Test1/Test1.summstats.txt (PREFORMATTED GWAS summary statistics)
# Test1/Test1.train.fam (training cohort sample IDs)
# Test1/Test1.train.phen (training cohort phenotypes)
# Test1/chrom1.Test1.train.dosage.gz (training cohort genotypes)
# Test1/chrom2.Test1.train.dosage.gz (training cohort genotypes)
# ... 
# Test1/Test1.val.fam (validation cohort sample IDs)
# Test1/Test1.val.phen (validation cohort phenotypes)
# Test1/chrom1.Test1.val.dosage.gz (validation cohort genotypes)
# Test1/chrom2.Test1.val.dosage.gz (validation cohort genotypes)
# ... 
```

## Running NPS
Test set #1 is small enough to run on desktop computers (MacOS and Linux are supported) without parallel processing. We provide a wrapper script (`run_all_chroms.sh`) to drive cluster job scripts sequentially on a desktop, by processing one chromosome at a time. To run test set #1 on computer clusters, see the [SGE instruction](https://github.com/sgchun/nps/blob/master/SGE.md) and [LSF instruction](https://github.com/sgchun/nps/blob/master/LSF.md). 

1. **Standardize genotypes.** This step standardizes the genotype data of training cohort to the mean of 0 and variance of 1. The training genotype files should be in the [dosage format](https://github.com/sgchun/nps/blob/master/FileFormats.md) and named as chrom*N*.*DatasetID*.dosage.gz. The command arguments are:
    * (1) directory where training genotype files are: `testdata/Test1`
    * (2) *DatasetID* of training genotype files: `Test1.train`. 
   ```bash
   cd nps-1.1.1/
   
   ./run_all_chroms.sh sge/nps_stdgt.job testdata/Test1 Test1.train
   ```

2. **Configure an NPS run.** For test set #1, which has ~100,000 genomewide SNPs, we recommend a window size of 80 SNPs. In general, for ~5,000,000 genome-wide SNPs we recommend to use 4,000-SNP windows. The command arguments are:
    * `--gwas`: GWAS summary statistics file (`testdata/Test1/Test1.summstats.txt`)
    * `--train-dir`: directory where training genotype files are (`testdata/Test1`)
    * `--train-dataset`: *DatasetID* of training cohort (`Test1.train`)
    * `--window-size`: genomic window size (`80`)
    * `--out`: directory to store NPS data (`testdata/Test1/npsdat`); all NPS output files are saved in this directory.
   ```bash
   Rscript npsR/nps_init.R --gwas testdata/Test1/Test1.summstats.txt \
       --train-dir testdata/Test1 \
       --train-dataset Test1.train \
       --window-size 80 \
       --out testdata/Test1/npsdat
   ```
   
   We provides the following optional arguments: 
    * `--train-fam`: sample file for training cohort in [FAM format](https://github.com/sgchun/nps/blob/master/FileFormats.md) (by default, {train dir}/{train dataset ID}.fam, e.g. `testdata/Test1/Test1.train.fam`)
    * `--train-phen`: phenotype file for training cohort. See [phenotype file format](https://github.com/sgchun/nps/blob/master/FileFormats.md). If this is not specified, NPS searches phenotype data first in the provided fam file. If phenotypes are missing in the fam file, NPS looks up a default phenotype file ({train dir}/{train dataset ID}.phen, e.g. `testdata/Test1/Test1.train.fam`).  

3. **Set up a special partition for GWAS-significant SNPs.** The command argument is: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
   ```bash
   ./run_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/
   ```

4. **Set up the decorrelated "eigenlocus" space.** This step projects genetic data into an orthogonalized domain and is the most time-consuming step in NPS. For the best performance, we recommend to run NPS on shifted overlapping windows and then to merge the results in the step (7). The recommended window shifts are 0, *WindowSize* \* 1/4, *WindowSize* \* 2/4 and WindowSize \* 3/4 SNPs, where *WindowSize* is the size of genomic windows. For test set #1, set the *WindowSize* to 80, thus the window shifts should be `0`, `20`, `40` and `60`. For the default window size of `4000`, we recommend the window shifts of `0`, `1000`, `2000` and `3000`. The command arguments are: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) window shift: `0`, `20`, `40` or `60`  
   ```bash
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 0
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 20
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 40
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 60
    ```
   
5. **Partition the rest of genome.** First, we define the partition cut-offs by running `nps_prep_part.R`. We recommend 10-by-10 double-partitioning on the intervals of eigenvalues of projection and estimated effect sizes in the eigenlocus space. The command arguments are:
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) Number of partitions on eigenvalues: `10`
    * (3) Number of partitions on estimated effects: `10` 
   ```
   Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 10 10 
   ```
   
   Then, we calculate partitioned polygenic risk scores in the training cohort by running `nps_part.job`. The command arguments are: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) window shift: `0`, `20`, `40` or `60`  
   ```
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 0
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 20
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 40
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 60
   ```

6. **Estimate shrinkage weights for each partition.** The command argument is: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
   If the sex is included in the FAM file of training data, NPS will automatically train the prediction model with the sex covariate.
   ```bash
   Rscript npsR/nps_reweight.R testdata/Test1/npsdat/ 
   ```
   
7. **Evaluate the accuracy of trained prediction model in a validation cohort.** First, NPS calculates polygenic risk scores for each individual in the validation cohort chromosome by chromosome. Use `nps_score.dosage.job` if genotype files are prepared in the dosage format (test set #1 and #2). If the validation cohort is prepared in Oxford bgen format, use `nps_score.bgen.job` instead. To use nps_score.bgen.job, [qctool version 2](https://www.well.ox.ac.uk/~gav/qctool_v2/) is required. nps_score.bgen.job calculates polygenic scores by running `qctool -risk-score`, and this requires that the SNP IDs are identical between training and validation genotype files. If needed, the SNP IDs of bgen files can be updated by running the `qctool -map-id-data` command. 

The command arguments for nps_score.dosage.job and nps_score.bgen.job are:
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) directory where validation cohort genotype files are: `testdata/Test1`
    * (3) *DatasetID* for validation genotypes files: `Test1.val`
    * (4) window shift: `0`, `20`, `40` or `60`
   NPS expects genotype file names of chrom*N*.*DatasetID*.dosage.gz for `nps_score.dosage.job` and chrom*N*.*DatasetID*.bgen for `nps_score.bgen.job`. 
   ```bash
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0 
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 20
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 40
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 60
   ```
   
   Then, NPS combines polygenic scores across all chromosomes and shifted windows with `nps_val.R` and reports the overall prediction accuracies. The command arguments are: 
   * `--out`: NPS data directory (`testdata/Test1/npsdat`)
   * `--val-dataset`: *DatasetID* for validation cohort: `Test1.val`
   ```  
   Rscript npsR/nps_val.R --out testdata/Test1/npsdat --val-dataset Test1.val 
   ```
   
   We provides the following optional arguments: 
    * `--val-dir`: directory of validation cohort data if it is different from training data dir
    * `--val-fam`: sample file for the validation cohort (by default, {val dir}/{val dataset ID}.fam, e.g. `testdata/Test1/Test1.val.fam`)
    * `--val-phen`: phenotype file for training cohort. If this is not specified, NPS searches phenotype data first in the provided fam file. If phenotypes are also missing in the fam file, NPS look up a default phenotype file ({val dir}/{val dataset ID}.phen, e.g. `testdata/Test1/Test1.val.fam`).
   See [File Formats](https://github.com/sgchun/nps/blob/master/FileFormats.md) for the details. If the risk prediction model was trained with the sex covariate at the step (6), NPS will incorporate the sex in the validation model as well using the sex in validation sample fam file.
   
   NPS stores the polygenic risk scores computed for the validation cohort in `testdata/Test1/Test1.val.phen.nps_score`. For test set #1, NPS will print out the following accuracy statistics: 
   > ...  
   > Polygenic scores are saved in testdata/Test1/Test1.val.phen.nps_score.  
   > ...   
   > Area under the curve: 0.8776  
   > 95% CI: 0.8589-0.8963 (DeLong)  
   > ...  
   > Nagelkerke's R2 = 0.3176188  
   > Tail OR (5%): 15.56794  
   > Done  

### Running NPS on test set #2
The steps to run NPS on test set #2 is simlar to that of test set #1. There are only a few small differences: With ~5 million genomewide SNPs, the genomic window size should be set to 4,000 (default value). And with the window size of 4,000, we recommend to run NPS on overlapping windows shifted by 0, 1,000, 2,000, and 3,000 SNPs. For step-by-step instructions, see [SGE](https://github.com/sgchun/nps/blob/master/SGE.md) and [LSF](https://github.com/sgchun/nps/blob/master/LSF.md) sections. We do not recommend running NPS on genomewide datasets using desktop computers.

NPS reports the following prediction accuracy with test set #2: 
> ...  
> Polygenic scores are saved in testdata/Test2/Test2.val.phen.nps_score.  
> ...  
> Area under the curve: 0.7964  
> 95% CI: 0.7701-0.8227 (DeLong)  
> ...  
> Nagelkerke's R2 = 0.1736275  
> Tail OR (5%): 7.163156  
> Done  

## Citation
> Chun et al. Non-parametric polygenic risk prediction using partitioned GWAS summary statistics.  
> BioRxiv 2020. doi: 10.1101/370064 (preprint).
