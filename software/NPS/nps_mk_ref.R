#! /usr/bin/env Rscript
rm(list=ls())

library(bigreadr)
library(optparse)

# Parameter setting
args_list = list(
  make_option("--val", type="character", default=NULL,
              help="INPUT: val data", metavar="character"),
  make_option("--valpheno", type="character", default=NULL,
              help="INPUT: valpheno data", metavar="character"),
  make_option("--outpath", type="character", default=NULL,
              help="OUTPUT: outpath", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

fam <- fread2(paste0(opt$val, ".fam"))
pheno <- fread2(opt$valpheno)[, 1]
pheno <- ifelse(is.na(pheno), -9, pheno)
fam_p <- cbind(fam[, -6], pheno)
write.table(fam_p, file = paste0(opt$outpath, "ref.fam"), 
            col.names = F, row.names = F, quote = F)
