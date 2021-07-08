## Input files for NPS
To run NPS, you need the following set of files: 

1. **GWAS summary statistics.** This is a *tab-delimited* text file format with the following seven required columns: 
     - **chr**: chromosome name starting with "chr." Currently, NPS expects only chromosomes 1-22. Chromosome names should be designated by "chr1", ..., "chr22".
     - **pos**: base position of SNP.
     - **ref** and **alt**: reference and alternative alleles of SNP, respectively.
     - **reffreq**: allele frequency of reference allele in the discovery GWAS cohort. 
     - **pval**: p-value of association. 
     - **effalt**: estimated *per-allele* effect size of *the alternative allele*. For case/control GWAS, log(Odds Ratio) should be used. NPS will convert **effalt** to effect sizes relative to *the standardized genotype* internally using **reffreq**.  
     ```
     chr	pos	ref	alt	reffreq	pval	effalt
     chr1	11008	C	G	0.9041	0.1126	-0.0251
     chr1	11012	C	G	0.9041	0.1126	-0.0251
     chr1	13116	T	G	0.8307	0.615	0.0071
     chr1	13118	A	G	0.8307	0.615	0.0071
     chr1	14464	A	T	0.8386	0.476	-0.0105
     ...
     ```

2. **Training genotype files.** Training genotype files should be in the qctool dosage format and named as "chrom*N*.*DatasetID*.dosage.gz" for each chromosome. Genotype files in bgen format can be converted to the dosage files by running [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/examples/converting.html) with the `-ofiletype dosage` option. NPS allows only *biallelic* variants. 
   
3. **Training sample file (.fam).** Sample information of training cohort should be provided in [PLINK FAM format](https://www.cog-genomics.org/plink2/formats#fam). The samples in the .fam file should appear in the exactly same order as in the genotype dosage files. The sex of sample (5-th column) is optional ("0" or "-9" for missing; "1" for male; "2" for female). If the sex is provided, NPS will incorporate the sex covariate in the PRS model. The 6-th column is for phenotype data and can be specified here or in a separeate phenotype file. 
   ```
   trainF2  trainI2  0  0  1 -9
   trainF3  trainI3  0  0  2 -9
   trainF39 trainI39 0  0  1 -9
   trainF41 trainI41 0  0  2 -9
   trainF58 trainI58 0  0  1 -9
   ```
4. **Training phenotype file (.phen).** Phenotypes of the .fam file can be overridden by a .phen file (use `nps_init.R --train-phen` option). This is a tab-delimited file with three columns: "FID", "IID", and "Outcome". FID and IID correspond to the family and individual IDs in the .fam file. The name of phenotype should be "Outcome". Binary phenotypes (case/control) are specified by "1" and "2", respectively; "0" and "-9" denote missing phenotype. For quantitative phenotypes, "-9" represents a missing phenotype value. 
   ```
   FID   IID    Outcome
   trainF2  trainI2  1
   trainF39 trainI39 1
   trainF3  trainI3  2
   trainF41 trainI41 2
   trainF58 trainI58 1
   ```
5. **Validation genotype files.** Validation genotypes can be in the dosage or .bgen format. If they are in .bgen format, the files should be named as "chrom*N*.*DatasetID*.bgen". 

6. **Validation sample file (.fam).** Similar to the training .fam file. 

7. **Training phenotype file (.phen).** Similar to training .phen file. 

