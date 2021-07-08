#! /usr/bin/env Rscript

rm(list=ls())
library(bigsnpr)
library(optparse)

## Input parameters
args_list = list(
  make_option("--mis", type="character", default=NULL,
              help="INPUT: missing data", metavar="character"), 
  make_option("--imp", type="numeric", default=NULL,
              help="INPUT: imputation data", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(mis = "/net/mulan/disk2/yasheng/predictionProject/plink_file/hm3/chr22", 
#             imp = "/net/mulan/disk2/yasheng/predictionProject/plink_file/hm3/chr22_imp")

## fast imputation
if(file.exists(paste0(opt$mis, ".bk")) == F){
  mis <- snp_readBed(paste0(opt$mis, ".bed"))
}
mis_bed <- snp_attach(paste0(opt$mis, ".rds"))
mis_bed$genotypes <- snp_fastImputeSimple(mis_bed$genotypes)
snp_writeBed(mis_bed, paste0(opt$imp, ".bed"))

## 
system(paste0("rm ", opt$mis, ".rds"))
system(paste0("rm ", opt$mis, ".bk"))
