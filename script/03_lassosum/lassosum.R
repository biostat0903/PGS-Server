#! /usr/bin/env Rscript
library(plyr)
library(bigreadr)
library(optparse)
library(lassosum)
library(doParallel)

## Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data prefix", metavar="character"), 
  make_option("--valid_genotype", type="character", default=NULL,
              help="INPUT: validation genotype prefix", 
              metavar="character"),
  make_option("--valid_phenotype", type="character", default=NULL,
              help="INPUT: validation phenotype prefix", 
              metavar="character"),
  make_option("--n", type="numeric", default=NULL,
              help="INPUT: sample size", metavar="character"),
  make_option("--population", type="character", default=NULL,
              help="INPUT: population", metavar="character"), 
  make_option("--covariates", type="character", default=NULL,
              help="INPUT: covariates", metavar="character"), 
  make_option("--outpath", type="character", default=NULL,
              help="OUTPUT: outpath", metavar="character"),
  make_option("--thread", type="numeric", default=NULL,
              help="OUTPUT: thread", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## validation data
bim_file <- fread2(paste0(opt$valid_genotype, ".bim"))
fam_file <- fread2(paste0(opt$valid_genotype, ".fam"))
pheno_file <- fread2(opt$valid_phenotype)[, 1]

if (length(unique(bim_file[, 1])) != 22){
  stop("ERROR: We need validation genotype from chr1-chr22!")
}
if (nrow(fam_file) != length(pheno_file)){
  stop("ERROR: Validation phenotype is not matched to validation genotype!")
}
if(is.null(opt$covariate) == F){
  cov_file <- fread2(opt$covariates)
  if (nrow(fam_file) != nrow(cov_file)){
    stop("ERROR: Validation phenotype is not matched to validation covariates!")
  }
}
idx_val <- fam_file[, 1]

## summary statistics input
summstats <- fread2(opt$summ, select =  c(1, 2, 3, 6, 7, 9, 10))
colnames(summstats) <- c("chr", "snp", "pos", "a1", "a2", "beta", "se")

## p to cor
t <- summstats$beta/summstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
p_val_z <- ifelse(p_val == 0, min(p_val[-which(p_val==0)]), p_val)
cor <- p2cor(p = p_val_z, n = opt$n, sign = summstats$beta)

## block information
setwd(system.file("data", package="lassosum"))
LDblocks <- paste0(opt$population, ".hg19")

# parallel setting
threads <- as.numeric(opt$thread)
cl <- makeCluster(threads, type="FORK")

## estimation
out <- lassosum.pipeline(cor = cor, chr = summstats$chr, pos = summstats$pos,
                         A1 = summstats$a1, A2 = summstats$a2, destandardize = F,
                         ref.bfile = opt$valid_genotype, test.bfile = opt$valid_genotype,
                         LDblocks = LDblocks,  cluster = cl)

## validation
if (is.null(opt$covariates)){
  pheno_val <- data.frame(FID = idx_val, IID = idx_val, pheno_file)
  valid <- validate(out, pheno = pheno_val, plot = F)
} else {
  pheno_val <- data.frame(FID = idx_val, IID = idx_val, pheno_file)
  valid <- validate(out, pheno = pheno_val, covar = cov_file, plot = F)
}

## best
out_best <- subset(out, s = valid$best.s, lambda = valid$best.lambda)
esteff_best <- data.frame(snp = summstats$snp[out_best$sumstats$order], 
                          a1 = summstats$a1[out_best$sumstats$order], 
                          beta = out_best$beta)
esteff_best_nz <- esteff_best[esteff_best[, 3] != 0, ]

write.table(esteff_best_nz, file = paste0(opt$outpath, "lassosum_esteff.txt"), 
            row.names = F, col.names = F, quote = F)
