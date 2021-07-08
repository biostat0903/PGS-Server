#! /usr/bin/env Rscript
rm(list=ls())
library(plyr)
library(bigreadr)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data", metavar="character"),
  make_option("--valpath", type="character", default=NULL,
              help="INPUT: val path", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/home/yasheng/comprsWeb/example_data/all/summary", 
#             valid = "/home/yasheng/comprsWeb/example_data/val/valid")

summstats <- fread2(paste0(opt$summ, ".assoc.txt"), 
                    select = c(2, 3, 6, 7, 8, 9, 10, 11))
colnames(summstats) <- c("rsid", "pos", "alt", "ref", "maf", "effalt", "se", "pval")
t <- summstats$effalt/summstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
p_val_z <- ifelse(p_val == 0, min(p_val[-which(p_val==0)]), p_val)
summstats$reffreq <- 1 - summstats$maf
summstats$pval <- p_val_z
summstats$maf <- NULL
summstats$se <- NULL

snpinfo <- fread2(paste0(opt$valpath, "merge.ref.snpinfo"))
snpinfo <- snpinfo[snpinfo[, 1] != "chromosome", ]
summstats_merge <- merge(snpinfo, summstats, by = "rsid", all.x = T)
if (any(is.na(summstats_merge$ref))){
  summstats_merge[is.na(summstats_merge$ref), 10] <- 0
  summstats_merge[is.na(summstats_merge$ref), 11] <- 1
  summstats_merge[is.na(summstats_merge$ref), 12] <- 0.5
  summstats_merge[is.na(summstats_merge$ref), 9] <- summstats_merge$alleleA[is.na(summstats_merge$ref)]
}
summstats_merge$effalt <- ifelse(summstats_merge$ref == summstats_merge$alleleA, 
                                 summstats_merge$effalt, 
                                 -summstats_merge$effalt)
summstats_inter <- summstats_merge[, c(2, 4, 5, 6, 12, 11, 10)]
summstats_inter[, 1] <- paste0("chr", summstats_inter[, 1])
colnames(summstats_inter) <- c("chr", "pos", "ref", "alt", "reffreq", "pval", "effalt")
summstats_inter <- summstats_inter[order(summstats_inter$chr), ]

write.table(summstats_inter, file = paste0(opt$summ, ".nps"), 
            row.names = F, quote = F, sep = "\t")
