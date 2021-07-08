rm(list=ls())
library(bigreadr)
library(optparse)

# Parameter setting
args_list = list(
  make_option("--valpath", type="character", default=NULL,
              help="INPUT: val data", metavar="character"),
  make_option("--windowsize", type="numeric", default=NULL,
              help="INPUT: valpheno data", metavar="character"),
  make_option("--tmppath", type="character", default=NULL,
              help="INPUT: valpheno data", metavar="character"),
  make_option("--outpath", type="character", default=NULL,
              help="OUTPUT: outpath", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

win_str <- as.integer(opt$windowsize*c(0, 0.25, 0.5, 0.75))
cat(win_str,  "\n")
for (chr in 1: 22){

  snp_info_str <- paste0(opt$valpath, "chrom", chr, ".ref.snpinfo")
  snp_info <- fread2(snp_info_str)

  tail_mat <- matrix(NA, nrow(snp_info), 4)
  for (winshift in 1: length(win_str)){
    tail_str <- paste0(opt$tmppath, "/ref.win_", win_str[winshift],
                       ".adjbetahat_tail.chrom", chr, ".txt")
    tail_mat[, winshift] <- fread2(tail_str)[, 1]
  }
  pg_mat <- matrix(NA, nrow(snp_info), 4)
  for (winshift in 1: length(win_str)){
    pg_str <- paste0(opt$tmppath, "/ref.win_", win_str[winshift],
                     ".adjbetahat_pg.chrom", chr, ".txt")
    pg_mat[, winshift] <- fread2(tail_str)[, 1]
  }

  eff_mat <- rowSums(pg_mat) + rowSums(tail_mat)
  eff_dat <- data.frame(snp_info[, 2], snp_info[, 5], eff_mat)
  write.table(eff_dat, file = paste0(outpath, "nps_esteff", "_chr", chr, ".txt"),
              col.names = F, row.names = F, quote = F)
}

