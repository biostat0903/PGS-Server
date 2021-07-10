#! /usr/bin/env Rscript
library(plyr)
library(bigreadr)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--plink_str", type="character", default=NULL,
              help="INPUT: validation genotype directory",
              metavar="character"),
  make_option("--valid_genotype", type="character", default=NULL,
              help="INPUT: validation genotype directory",
              metavar="character"),
  make_option("--valid_phenotype", type="character", default=NULL,
              help="INPUT: validation phenotype directory",
              metavar="character"),
  make_option("--esteff_path", type="character", default=NULL,
              help="INPUT: estimated effect size directory",
              metavar="character"),
  make_option("--covariates", type="character", default=NULL,
              help="INPUT: covariates file", metavar="character"),
  make_option("--index", type="character", default=NULL,
              help="INPUT: index", metavar="character"),
  make_option("--outpath", type="character", default=NULL,
              help="INPUT: outpath", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(plink_str = "/public/home/biostat03/project/compProject/COM-PGS-main/software/plink",
#             valid_genotype = "/public/home/biostat03/project/compProject/COM-PGS-main/example_data/val/valid", 
#             valid_phenotype = "/public/home/biostat03/project/compProject/COM-PGS-main/example_data/val/valid_pheno.txt", 
#             esteff_path = "/public/home/biostat03/project/compProject/COM-PGS-main/example_data/output/PRSCS/", 
#             index = "r2", 
#             outpath = "/public/home/biostat03/project/compProject/COM-PGS-main/example_data/output/PRSCS/")

pheno <- fread2(opt$valid_phenotype)[, 1]
if (!is.null(opt$covariates)){
  cov <- as.matrix(fread2(opt$covariates))
  cov_dat <- data.frame(y = pheno, cov)
  cov_eff <- predict(lm(y~., data = cov_dat))
} 
# phi_str <- c("1e-06", "1e-04", "1e-02", "1e+00")
phi_str <- c("1e-02", "1e+00")
res_mat <- vector("numeric", length(phi_str))

## select phi
for (phi in 1: length(phi_str)){
  val_pheno <- matrix(NA, length(pheno), 22)
  for (chr in 21: 22){
    esteff_str <- paste0(opt$outpath, "/esteff_phi", phi_str[phi], "_chr", chr, ".txt")
    if(file.exists(esteff_str)){
      val_str <- paste0(opt$outpath,"/val_phi", phi_str[phi], "_chr", chr)
      val_cmd <- paste0(opt$plink_str, " --silent --bfile ", opt$valid_genotype, " --score ", esteff_str, " 2 4 6 sum ", 
                        " --out ", val_str)
      system(val_cmd)
      system(paste0("rm ", val_str, ".log"))
      val_pheno_chr <- fread2(paste0(val_str, ".profile"))
      val_pheno[, chr] <- val_pheno_chr[, 6]
    } else {
      val_pheno[, chr] <- 0
    }
  }
  val_pheno_sum <- rowSums(val_pheno, na.rm = T) 
  if (opt$index == "auc"){
    val_pheno_cov <- val_pheno_sum[!is.na(pheno)] + cov_eff[!is.na(pheno)]
    res_mat[phi] <- pROC::roc(pheno[!is.na(pheno)], val_pheno_cov)$auc
  } else {
    res_mat[phi] <- cor(pheno[!is.na(pheno)], val_pheno_sum[!is.na(pheno)])^2
  }
} 
best_phi <- phi_str[which.max(phi_str)]
cat(best_phi, "\n")

for (chr in 21: 22){
  
  ## effect file
  bestphi_est_str <- paste0(opt$outpath, "/esteff_bestphi",
                            "_chr", chr, ".txt")
  mv_est_cmd <- paste0("mv ", opt$outpath, "/esteff_phi", best_phi,
                       "_chr", chr, ".txt ", bestphi_est_str)
  system(mv_est_cmd)
  system(paste0("gzip ", bestphi_est_str))
  
  ## remove other files
  system(paste0("rm ", opt$outpath, "/esteff_phi*",
                "_chr", chr, ".txt "))
  
  ## validation file
  bestphi_val_str <- paste0(opt$outpath, "/val_bestphi", 
                            "_chr", chr, ".profile")
  mv_val_cmd <- paste0("mv ", opt$outpath, "/val_phi", best_phi,
                       "_chr", chr, ".profile ", 
                       bestphi_val_str)
  system(mv_val_cmd)
  system(paste0("gzip ", bestphi_val_str))
  
  ## remove other validation files
  system(paste0("rm ", opt$outpath, "/val_phi*",
                "_chr", chr, ".profile"))
}

write.table(res_mat, file = paste0(opt$outpath, "/", opt$index, ".txt"), 
            row.names = F, col.names = F, quote = F)
