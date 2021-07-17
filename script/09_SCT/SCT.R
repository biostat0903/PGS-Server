#! /usr/bin/env Rscript
rm(list=ls())
library(plyr)
library(bigsnpr)
library(bigreadr)
library(tidyverse)
library(optparse)

## Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file directory",
              metavar="character"),
  make_option("--valid_genotype", type="character", default=NULL,
              help="INPUT: validation genotype directory",
              metavar="character"),
  make_option("--valid_phenotype", type="character", default=NULL,
              help="INPUT: validation phenotype directory",
              metavar="character"),
  make_option("--p_len", type="numeric", default=10,
              help="INPUT: P value length", metavar="character"),
  make_option("--r2_val", type="character", default=c(0.01, 0.1, 0.2, 0.5, 0.8),
              help="INPUT: r2 values", metavar="character"),
  make_option("--dist_str", type="character", default=c(50, 100, 200, 500),
              help="INPUT: window sizes", metavar="character"),
  make_option("--covarites", type="character", default="cov",
              help="INPUT: covariates", metavar="character"),
  make_option("--thread", type="numeric", default=1,
              help="INPUT: thread number", metavar="character"), 
  make_option("--outpath", type="character", default=NULL,
              help="OUTPUT: output path", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

### Parameters
dist <- unlist(strsplit(opt$dist_str, ","))
r2 <- unlist(strsplit(opt$r2_val, ","))

### validation data
bim_file <- fread2(paste0(opt$valid_genotype, ".bim"))
fam_file <- fread2(paste0(opt$valid_genotype, ".fam"))
if (length(unique(bim_file[,1])) != 22){
  stop("We need validation genotype from chr1-chr22!")
}
if(file.exists(paste0(opt$valid_genotype, ".bk")) == F){
  val_bed <- snp_readBed(paste0(opt$valid_genotype, ".bed"))
}
val_bed <- snp_attach(paste0(opt$valid_genotype, ".rds"))
y <- fread2(opt$valid_phenotype)[, 1]
G <- snp_fastImputeSimple(val_bed$genotypes)
CHR <- val_bed$map$chromosome
POS <- val_bed$map$physical.pos
n_snp <- dim(G)[2]
if (dim(G)[1] != nrow(fam_file)) {
  stop("ERROR: Validation phenotype is not matched to validation genotype!")
}

## summary statistics
summstats <- fread2(opt$summ, select =  c(1, 2, 3, 7, 6, 9, 10))
colnames(summstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "se") # calculate P
if (length(unique(summstats$chr)) != 22){
  stop("ERROR: We need summary statistics from chr1-chr22!")
}
t <- summstats$beta/summstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
summstats$pval <- ifelse(p_val == 0,
                         min(p_val[-which(p_val==0)]),
                         p_val)

### Map summary statistics and reference panel
map <- val_bed$map[, -3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
info_snp <- snp_match(summstats, map, strand_flip = F)
beta <- rep(0, n_snp)
lp_val <- rep(0, n_snp)
beta[map[, 2]%in%info_snp[, 5]] <- info_snp$beta
lp_val[map[, 2]%in%info_snp[, 5]] <- -log10(info_snp$pval)

## clumping
all_keep <- snp_grid_clumping(G,
                              infos.chr = CHR,
                              infos.pos = POS,
                              lpS = lp_val,
                              grid.thr.r2 = as.numeric(r2),
                              grid.base.size = as.numeric(dist),
                              ncores = opt$thread)

## threshold
bk_num <- as.numeric(as.POSIXct(Sys.time()))
bk_str <- paste0(opt$outpath, "SCT_est", bk_num)
multi_PRS <- snp_grid_PRS(G,
                          all_keep,
                          betas = beta,
                          lpS = lp_val,
                          backingfile = bk_str,
                          n_thr_lpS = opt$p_len,
                          ncores = opt$thread)
nn <- nrow(attr(all_keep, "grid"))
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = c(1:nn)) %>%
  unnest(cols = "thr.lp")
s <- nrow(grid2)

## subsample phenotype
if (is.null(opt$covariates)){

  y <- y[!is.na(y)]
  final_mod <- snp_grid_stacking(multi_PRS,
                                 y,
                                 ncores = opt$thread,
                                 K = 10)
} else{
  y <- fread2(opt$valid_phenotype)[, 1]
  covar <- fread2(opt$valid_phenotype)[!is.na(y), ]
  covar <- covar[!is.na(y)]
  y <- y[!is.na(y)]
  final_mod <- snp_grid_stacking(multi_PRS,
                                 y,
                                 covar.train = as.matrix(covar),
                                 ncores = opt$thread,
                                 K = 10)
}

## SCT
new_beta <- final_mod$beta.G
idx <- which(new_beta != 0)
snp_sig_SCT <- data.frame(map$rsid[idx],
                          map$a1[idx],
                          new_beta[idx])

# output
write.table(snp_sig_SCT, file = paste0(opt$outpath, "SCT_esteff.txt"),
            col.names = F, row.names = F, quote = F)
system(paste0("rm ", opt$valid_genotype, ".bk"))
system(paste0("rm ", bk_str, ".bk"))
system(paste0("rm ", bk_str, ".rds"))
