#! /usr/bin/env Rscript
rm(list=ls())
library(plyr)
library(bigsnpr)
library(bigreadr)
library(optparse)
library(tidyverse)

## Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--valid_genotype", type="character", default=NULL,
              help="INPUT: validation genotype directory",
              metavar="character"),
  make_option("--valid_phenotype", type="character", default=NULL,
              help="INPUT: validation phenotype directory",
              metavar="character"),
  make_option("--dist", type = "numeric", default = 500,
              help = "INPUT: the distance of LD matrix", 
              metavar = "character"),
  make_option("--covariates", type="character", default="cov",
              help="INPUT: covariates", metavar="character"),
  make_option("--chr", type="numeric", default=NULL,
              help="INPUT: chr number", metavar="character"),
  make_option("--outpath", type="character", default=NULL,
              help="OUTPUT: CT path", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## parameter
p_len <- 21

### validation data
bim_file <- fread2(paste0(opt$valid_genotype, ".bim"))
fam_file <- fread2(paste0(opt$valid_genotype, ".fam"))
if(file.exists(paste0(opt$valid_genotype, ".bk")) == F){
  val_bed <- snp_readBed(paste0(opt$valid_genotype, ".bed"))
}
val_bed <- snp_attach(paste0(opt$valid_genotype, ".rds"))

### covariates
if (opt$covariates == "cov"){
 
  y <- fread2(opt$valid_phenotype)[, 1]
} else {
  
  y <- fread2(opt$valid_phenotype)[, 1]
  covar <- fread2(opt$cov)[!is.na(y), ]
  covar <- covar[!is.na(y)]
  y <- y[!is.na(y)]
}

if (all(is.na(y)) == F){
  idx <- which(!is.na(y))
  sub_bk_str <- paste0(opt$valid_genotype, "_sub-",
                       as.numeric(as.POSIXlt(Sys.time())))
  val_bed <- snp_attach(snp_subset(val_bed, ind.row = idx, backingfile = sub_bk_str))
  y <- y[idx]
}

### genotype
G <- snp_fastImputeSimple(val_bed$genotypes)
CHR <- val_bed$map$chromosome
POS <- val_bed$map$physical.pos
n_snp <- dim(G)[2]
map <- val_bed$map[-(2:3)]
names(map) <- c("chr", "pos", "a1", "a0")

## summary data
summstats <- fread2(opt$summ,
                    select =  c(1, 2, 3, 5, 6, 7, 8, 9, 10))
colnames(summstats) <- c("chr", "rsid", "pos", "n_obs", 
                         "a1", "a0", "MAF", "beta", "beta_se")
summstats$n_eff <- summstats$n_obs
summstats$n_obs <- NULL
summstats$sgenosd <- 2*summstats$MAF*(1-summstats$MAF)
summstats$MAF <- NULL
info_snp <- snp_match(summstats, map, match.min.prop = 0.05)

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
info_chr <- info_snp$`_NUM_ID_`

## LDpred2_inf
corr <- snp_cor(G, ind.col = info_chr, 
                size = opt$dist)
ldsc <- snp_ldsc2(corr, df_beta)
h2_est <- ldsc[["h2"]]
cat (opt$chr, ":", h2_est, "\n")
if(h2_est < 0){
  beta_LDpred2 <- NA
} else {
  corr_sp <- bigsparser::as_SFBM(as(corr, "dgCMatrix"))
  beta_inf <- snp_ldpred2_inf(corr_sp, df_beta, h2 = h2_est)
  cat ("Inf model is ok!\n")
  
  ## LDpred2-auto
  auto <- snp_ldpred2_auto(corr_sp, df_beta, h2_init = h2_est)
  beta_auto <- auto[[1]]$beta_est
  beta_auto <- ifelse(is.na(beta_auto), 0, beta_auto)
  cat ("LDpred-auto model is ok!\n")
  
  ## LDpred2
  p_seq <- signif(seq_log(1e-4, 1, length.out = p_len), 2)
  h_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
  params <- expand.grid(p = p_seq, h2 = h_seq, sparse = c(FALSE, TRUE))
  beta_grid <- snp_ldpred2_grid(corr_sp, df_beta, params)
  idx_sp_na <- is.na(beta_grid[1, c(1: p_len)]) | abs(beta_grid[1, c(1: p_len)])>1
  idx_nosp_na <- is.na(beta_grid[1, (p_len+c(1: p_len))]) | abs(beta_grid[1, (p_len+c(1: p_len))])>1
  
  if (all(idx_nosp_na) & all(idx_sp_na)){
    ## output
    beta_LDpred2 <- data.frame(info_snp$rsid,
                               info_snp$a1,
                               beta_inf,
                               0,
                               0,
                               beta_auto)
  } else {
    if(sum(idx_sp_na) >= sum(idx_nosp_na)){
      idx_na <- c(which(idx_sp_na==F), (p_len + which(idx_sp_na==F)))
    } else {
      idx_na <- c(which(idx_nosp_na==F), (p_len + which(idx_nosp_na==F)))
    }
    
    beta_grid_na <- beta_grid[, idx_na]
    params_na <- params[idx_na, ]
    pred_grid <- big_prodMat(G, beta_grid_na, ind.col = info_chr)
    
    if (opt$covariates == "cov"){
      params_na[c("coef", "score")] <-
        big_univLinReg(as_FBM(pred_grid), y)[c("estim", "score")]
      params_na$idx_val  <- apply(pred_grid, 2, function(a) cor(a, y)^2)
    } else {
      params_na[c("coef", "score")] <-
        big_univLinReg(as_FBM(pred_grid), 
                       y, 
                       covar.train = as.matrix(covar))[c("estim", "score")]
      params_na$idx_val  <- apply(pred_grid, 2, AUC, target = y)
    }
    
    ## parameters
    params_na %>%
      mutate(sparsity = colMeans(beta_grid_na == 0), id = c(1: nrow(params_na))) %>%
      arrange(desc(score)) %>%
      mutate_at(4:8, signif, digits = 3)
    # no-sparsity effect
    best_grid_nosp <- params_na %>%
      mutate(id = c(1: nrow(params_na))) %>%
      filter(!sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      { beta_grid_na[, .$id] * .$coef }
    # sparsity effect
    best_grid_sp <- params_na %>%
      mutate(id = c(1: nrow(params_na))) %>%
      filter(sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      { beta_grid_na[, .$id] * .$coef }
    cat ("LDpred model is ok!\n")
    ## output
    beta_LDpred2 <- data.frame(info_snp$rsid,
                               info_snp$a1,
                               beta_inf,
                               best_grid_nosp,
                               best_grid_sp,
                               beta_auto)
  }
}

system(paste0("rm ", sub_bk_str, ".bk"))
system(paste0("rm ", sub_bk_str, ".rds"))

write.table(beta_LDpred2, 
            file = paste0(opt$outpath, "LDpred2_esteff_chr", opt$chr, ".txt"),
            col.names = F, row.names = F, quote = F)