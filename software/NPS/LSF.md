## Running NPS on LSF

### Before you run NPS
The cluster job scripts are located under [nps-1.1.1/sge/](https://github.com/sgchun/nps/tree/master/sge). These scripts run not only with SGE but also with UGER, LSF and Slurm schedulers. Depending on the system, it may be necessary to modify the provided job scripts to load required modules, for example:
```bash
###
# ADD CODES TO LOAD MODULES HERE
# ---------------------- EXAMPLE ----------------------------
# On clusters running environment modules and providing R-mkl
module add gcc/5.3.0 
module add R-mkl/3.3.2
# -----------------------------------------------------------
...
```

### nps_check command
We provide a command line tool (`nps_check.sh`) to verify the integrity of each job. We strongly recommend to check the completion of all jobs by running this tool before moving onto the next step. If nps_check detects an error, "FAIL" message will be printed. Otherwise, only "OK" messages will be reported. For example, running nps_check after nps_init reports the following checks:  
> NPS data directory: testdata/Test1/npsdat/  
> Verifying nps_init:  
> Checking testdata/Test1/npsdat//args.RDS ...OK (version 1.1)  
> Checking testdata/Test1/npsdat//log ...OK  
> Verifying nps_stdgt:  
> Checking testdata/Test1/chrom1.Test1.train ...OK  
> Checking testdata/Test1/chrom2.Test1.train ...OK  
> Checking testdata/Test1/chrom3.Test1.train ...OK  
> ...  

### Step-by-step instructions for test set #1
All steps have to run in the top-level NPS directory (nps-1.1.1/). The option `-J jobname[1-22]` will run NPS jobs over all 22 chromosomes in parallel. Currently, the number of chromosomes in the genome is fixed to 22 and not modifiable.
```
cd nps-1.1.1/

# Standardize genotypes
bsub -J stdgt[1-22] sge/nps_stdgt.job testdata/Test1 Test1.train

# Configure
Rscript npsR/nps_init.R --gwas testdata/Test1/Test1.summstats.txt \
       --train-dir testdata/Test1 \
       --train-dataset Test1.train \
       --window-size 80 \
       --out testdata/Test1/npsdat

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Separate GWAS-significant SNPs
bsub -J gwassig[1-22] sge/nps_gwassig.job testdata/Test1/npsdat/

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Set up the eigenlocus space 
bsub -J decor[1-22] sge/nps_decor_prune.job testdata/Test1/npsdat/ 0 
bsub -J decor[1-22] sge/nps_decor_prune.job testdata/Test1/npsdat/ 20 
bsub -J decor[1-22] sge/nps_decor_prune.job testdata/Test1/npsdat/ 40 
bsub -J decor[1-22] sge/nps_decor_prune.job testdata/Test1/npsdat/ 60 

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Partition the rest of genetic variations
Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 10 10

# Calculate partitioned risk scores in the training cohort
bsub -J part[1-22] sge/nps_part.job testdata/Test1/npsdat/ 0
bsub -J part[1-22] sge/nps_part.job testdata/Test1/npsdat/ 20
bsub -J part[1-22] sge/nps_part.job testdata/Test1/npsdat/ 40
bsub -J part[1-22] sge/nps_part.job testdata/Test1/npsdat/ 60

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Estimate per-partition shrinkage weights
Rscript npsR/nps_reweight.R testdata/Test1/npsdat/

# Calculate polygenic scores for each chromosome and for each individual in the validation cohort
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0 
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 20 
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 40 
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 60 

# Check the results 
./nps_check.sh testdata/Test1/npsdat/ 

# Calculate overall polygenic scores and report prediction accuracies
Rscript npsR/nps_val.R --out testdata/Test1/npsdat --val-dataset Test1.val 
```

### Step-by-step instructions for test set #2
In most cases, 4GB memory space per a task will be sufficient for running NPS jobs. On LSF clusters, the memory requirement can be specified by `bsub -R 'rusage[mem=4000]'`. 

```bash
cd nps-1.1.1/

# Standardize genotypes
bsub -J stdgt[1-22] sge/nps_stdgt.job testdata/Test2/ Test2.train

# Configure
# CAUTION: This step needs at least 4G memory
# You need to run this as a job or open an interative session with "-R 'rusage[mem=4000]'"
Rscript npsR/nps_init.R --gwas testdata/Test2/Test2.summstats.txt \
    --train-dir testdata/Test2 \
    --train-dataset Test2.train \
    --out testdata/Test2/npsdat

# Check the results
./nps_check.sh testdata/Test2/npsdat/ 

# Separate GWAS-significant SNPs
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] sge/nps_gwassig.job testdata/Test2/npsdat/

# Check the results
./nps_check.sh testdata/Test2/npsdat/ 

# Set up the eigenlocus space
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 0
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 1000
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 2000
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Partition the rest of genetic variations
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 10 10

# Calculate partitioned risk scores in the training cohort
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 0
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 1000
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 2000
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Estimate per-partition shrinkage weights
Rscript npsR/nps_reweight.R testdata/Test2/npsdat/

# Calculate polygenic scores for each chromosome and for each individual in the validation cohort
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0 
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 1000   
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 2000   
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 3000   

# Check the results 
./nps_check.sh testdata/Test2/npsdat/ 

# Calculate overall polygenic scores and report prediction accuracies
Rscript npsR/nps_val.R --out testdata/Test2/npsdat --val-dataset Test2.val
```

