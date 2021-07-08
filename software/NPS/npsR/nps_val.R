VERSION <- "1.1"

cat("Non-Parametric Shrinkage", VERSION, "\n")

print.usage <- function() {
    cat("Usage:\n")
    cat("Rscript nps_val.R <work dir> <val dataset ID> <val fam file> <val phen file> [<WINSHIFT> ...]")
    cat("\nRscript [options] [<WINSHIFT> ...]\n")
    cat("    --out <work dir>\n")
    # cat("    --val-dataset <val dataset ID>\n")
    # cat("    --val-dir <val dir>\n")
    # cat("          default value: <training dataset dir>\n")
    # cat("    --val-fam <val fam file>\n")
    # cat("          default value: <val dir>/<val dataset ID>.fam\n")
    # cat("    --val-phen <val phen file>\n")
    # cat("          default value: Use <val fam file> or\n")
    # cat("          <val dir>/<val dataset ID>.phen\n")
    cat("\n\n")
}

# ASSERT <- function(test) {
#     if (length(test) == 0) {
#         stop(paste("ASSERT fail for empty conditional:",
#                    deparse(substitute(test))))
#     }
# 
#     if (is.na(test)) {
#         stop(paste("ASSERT fail for missing value:",
#                    deparse(substitute(test))))
#     }
#     
#     if (!test) {
#         stop(paste("ASSERT fail:", deparse(substitute(test))))
#     }
# }
# 
# #########################################################################
# 
cargs <- commandArgs(trailingOnly=TRUE)

#     
# if (any(startsWith(cargs, "--"))) {
#     ## Parse long options
# 
#     valdir <- NULL
#     valfamfile <- NULL
#     valphenofile <- NULL
#     valtag <- NULL
#     tempprefix <- NULL
#     
#     arg.index <- 1
# 
#     while (arg.index <= length(cargs)) {
# 
#         if (!startsWith(cargs[arg.index], "--")) {
#             ## start of window shifts
# 
#             if (any(startsWith(cargs[arg.index:length(cargs)], "--"))) {
#                 cat("Invalid WINSHIFTs: ",
#                     paste(cargs[arg.index:length(cargs)], collapse=", "),
#                     "\n\n")
#                 print.usage()
#                 q()
#             }
#             
#             break
#         }
# 
#         opt <- substr(cargs[arg.index], 3, nchar(cargs[arg.index]))
# 
#         if (opt == "val-dir") {
#             arg.index <- arg.index + 1
# 
#             if (arg.index > length(cargs)) {
#                 stop("missing option value for ", paste("--", opt, sep=''))
#             }
# 
#             valdir <- cargs[arg.index]
# 
#         } else if (opt == "val-fam") {
#             arg.index <- arg.index + 1
# 
#             if (arg.index > length(cargs)) {
#                 stop("missing option value for ", paste("--", opt, sep=''))
#             }
# 
#             valfamfile <- cargs[arg.index]
# 
#         } else if (opt == "val-phen") {
#             arg.index <- arg.index + 1
# 
#             if (arg.index > length(cargs)) {
#                 stop("missing option value for ", paste("--", opt, sep=''))
#             }
# 
#             valphenofile <- cargs[arg.index]
# 
#         } else if (opt == "val-dataset") {
#             arg.index <- arg.index + 1
# 
#             if (arg.index > length(cargs)) {
#                 stop("missing option value for ", paste("--", opt, sep=''))
#             }
# 
#             valtag <- cargs[arg.index]
# 
#         } else if (opt == "out") {
#             arg.index <- arg.index + 1
# 
#             if (arg.index > length(cargs)) {
#                 stop("missing option value for ", paste("--", opt, sep=''))
#             }
# 
#             tempprefix <- paste(cargs[arg.index], "/", sep='')
# 
#         } else {
#             cat("Cannot recognize the option:", cargs[arg.index], "\n\n")
#             print.usage()
#             q()
#         }
#         
#         arg.index <- arg.index + 1
#     }
#     
#     if (arg.index < length(cargs)) {
#         ## window shifts specified
#         WINSHIFT.list <- as.numeric(cargs[arg.index:length(cargs)])
# 
#         if (any(is.na(WINSHIFT.list))) {
#             stop("Invalid window shift (non-numeric): ", 
#                  paste(cargs[arg.index:length(cargs)], collapse=", "))
#         }
#     } else {
#         WINSHIFT.list <- NULL
#     }
# 
#     ## Check required values
#     if (is.null(tempprefix)) {
#         stop("--out required")
#     }
# 
#     if (is.null(valtag)) {
#         stop("--val-dataset required")
#     }
#         
# } else {
#     
#     if (length(cargs) < 4) {
#         print.usage()
#         q()
#     }
#     
#     tempprefix <- paste(cargs[1], "/", sep='')
# 
#     valtag <- cargs[2]
#     valfamfile <- cargs[3]
#     valphenofile <- cargs[4]
# 
#     if (length(cargs) > 4) {
#         
#         WINSHIFT.list <- as.numeric(cargs[5:length(cargs)])
#         
#         if (any(is.na(WINSHIFT.list))) {
#             stop("Invalid window shift (non-numeric): ", 
#                  paste(cargs[5:length(cargs)], collapse=", "))
#         }
#         
#     } else {
#         WINSHIFT.list <- NULL
#     }
# }
tempprefix <- cargs[2]
valtag <- cargs[4] 
samplesize <- cargs[6] 
WINSHIFT.list <- NULL
args <- readRDS(paste(tempprefix, "/args.RDS", sep=''))

traintag <- args[["traintag"]]
traindir <- args[["traindir"]]
WINSZ <- args[["WINSZ"]]

if (is.null(WINSHIFT.list)) {
    cat("Detecting window shifts :")

    part.files <- list.files(tempprefix, pattern="*.part.RDS")
        
    WINSHIFT.list <-
        sapply(part.files,
               function (s) strsplit(s, ".", fixed=TRUE)[[1]][1],
               simplify=TRUE)

    WINSHIFT.list <-
        sapply(WINSHIFT.list,
               function (s) strsplit(s, "_", fixed=TRUE)[[1]][2],
           simplify=TRUE)

    WINSHIFT.list <- sort(as.numeric(WINSHIFT.list))

    if (length(WINSHIFT.list) == 0) {
        cat("ERROR\n")
        stop("No window shift found")
    }

    cat(paste(WINSHIFT.list, collapse=" "), "\n")
}

## check if winshifts are valid
# if (any(is.na(WINSHIFT.list)) || any(WINSHIFT.list < 0) ||
#     any(WINSHIFT.list >= WINSZ)) {
# 
#     stop("Invalid shift (window size =", WINSZ, "): ",
#          paste(WINSHIFT.list, collapse=", "))
# }

## Set default parameters
# if (is.null(valdir)) {
#     valdir <- traindir
# }
# 
# if (is.null(valfamfile)) {
#     valfamfile <- paste(valdir, "/", valtag, ".fam", sep='')
#     cat("Using the default validation fam file path:", valfamfile, "\n")
# }
# 
# if (!file.exists(valfamfile)) {
#     stop("File does not exists:", valfamfile)
# }
# 

#########################################################################
### validation

# phenotypes
# vlfam <- read.delim(valfamfile, sep=" ", header=FALSE,
#                     stringsAsFactors=FALSE)
# 
# if (ncol(vlfam) != 6) {
#     # re-try with tab delimination
# 
#     vlfam <- read.delim(valfamfile, sep="\t", header=FALSE,
#                         stringsAsFactors=FALSE)
# }

# if (ncol(vlfam) != 6) {    
#     stop(valfamfile, " does not have standard 6 columns (space or tab-delimited)")
# }
# 
# if (any(duplicated(paste(vlfam[, 1], vlfam[, 2], sep=":")))) {
#     stop("Duplicated FID IID combinations:", valfamfile)
# }

# vlphen <- NULL

# if (is.null(valphenofile)) {
#     if (all(vlfam[, 6] == 0) || all(vlfam[, 6] == -9)) {
#         cat("Phenotype data are empty in", valfamfile, ".\n")
# 
#         ## default phenotype file path
#         valphenofile <- paste(valdir, "/", valtag, ".phen", sep='')
# 
#     } else {
#         cat("Use phenotypes in", valfamfile, ".\n")
#         vlphen <- vlfam[, c(1, 2, 6)]
#         colnames(vlphen) <- c("FID", "IID", "Outcome")
#     }
# }
# 
# if (is.null(vlphen)) {
#     cat("Reading", valphenofile, "for phenotype data.\n")
# 
#     if (!file.exists(valphenofile)) {
#         stop("File does not exists:", valphenofile)
#     }
# 
#     vlphen <- read.delim(valphenofile, sep="\t", header=TRUE,
#                          stringsAsFactors=FALSE)
# }

# if (length(intersect(colnames(vlphen), c("FID", "IID", "Outcome"))) != 3) {
#     stop(valphenofile, " does not include standard columns: FID IID Outcome (tab-delimited")
# }
# 
# if (any(duplicated(paste(vlphen$FID, vlphen$IID, sep=":")))) {
#     stop("Duplicated FID IID combinations:", valphenofile)
# }
# 
# rownames(vlphen) <- paste(vlphen$FID, vlphen$IID, sep=":")

# No sample IDs match 
## if (length(intersect(rownames(vlphen), paste(vlfam[, 1], vlfam[, 2], sep=":")))
##    == 0) {
##    stop("IID/FID does not match between .fam and phenotype files")
## }

# samples in .fam but not phenotype file
## missing.entry <- setdiff(paste(vlfam[, 1], vlfam[, 2], sep=" "),
##                            paste(vlphen$FID, vlphen$IID, sep=" "))

## if (length(missing.entry) > 0) {
##    cat("FID IID\n")
##    cat(paste(missing.entry, collapse="\n"))
##    cat("\n")
## #    stop("The above samples declared in ", valfamfile,
## #         " are missing in the phenotype file: ",
## #         valphenofile)
## }

# vlphen <- vlphen[paste(vlfam[, 1], vlfam[, 2], sep=":"), ]
# vlphen$FID <- vlfam[, 1]
# vlphen$IID <- vlfam[, 2]

## vlphen$Outcome[is.na(vlphen$Outcome)] <- -9

## ASSERT(all(vlphen$FID == vlfam[, 1]))
## ASSERT(all(vlphen$IID == vlfam[, 2]))

# if (!is.numeric(vlphen$Outcome)) {
#     stop("phenotype values are not numeric")
# }

# cat("Validation cohort:\n")
# cat("Size of validation cohort:", nrow(vlfam), "\n")
# 
# binary.phen <- TRUE
# 
# if (length(unique(vlphen$Outcome[!is.na(vlphen$Outcome)])) > 4) {
# 
#     vlphen$Outcome[which(vlphen$Outcome == -9)] <- NA
# 
#     cat("# of missing phenotypes:", sum(is.na(vlphen$Outcome)), "\n")
#     cat("# of non-missing phenotypes:", sum(!is.na(vlphen$Outcome)), "\n")
#  
#     binary.phen <- FALSE
# 
# } else {
#     ## Binary phenotype
# 
#     if (all(vlphen$Outcome == 0 | vlphen$Outcome == 1 |
#              vlphen$Outcome == -9, na.rm=TRUE)) {
#         ## NPS v1.0 used 0/1/-9 encoding
#         outcome[which(vlphen$Outcome == -9)] <- NA 
# 
#     } else {
#         ## 1/2/0/-9 encoding
# #        ASSERT(any(vlphen$Outcome == 1))
# #        ASSERT(any(vlphen$Outcome == 2))
# #        ASSERT(all(vlphen$Outcome == 1 | vlphen$Outcome == 2 |
# #                   vlphen$Outcome == -9 | vlphen$Outcome == 0, na.rm=TRUE))
# 
#         if (!all(vlphen$Outcome == 1 | vlphen$Outcome == 2 |
#                  vlphen$Outcome == -9 | vlphen$Outcome == 0, na.rm=TRUE)) {
#             stop("Binary phenotype values can be only 1, 2, 0, or -9")
#         } 
# 
#         ## for backward compatibility
#         ## Outcome code 2 -> 1 (case)
#         ## Outcome code 1 -> 0 (control)
#         ## Outcome code 0 -> -9 (missing)
#         ## Outcome code -9 -> -9 (missing)
#         
#         ## recode to 0/1/-9
#         outcome <- rep(NA, nrow(vlphen))
#         outcome[which(vlphen$Outcome == 2)] <- 1
#         outcome[which(vlphen$Outcome == 1)] <- 0
#         outcome[which(vlphen$Outcome == 0)] <- NA
#         outcome[which(vlphen$Outcome == -9)] <- NA 
#         vlphen$Outcome <- outcome
#     }
#     
#     cat("# of missing phenotypes:", sum(is.na(vlphen$Outcome)), "\n")
#     cat("# of cases:", sum(vlphen$Outcome == 1, na.rm=TRUE), "\n")
#     cat("# of controls:", sum(vlphen$Outcome == 0, na.rm=TRUE), "\n")
# 
#     binary.phen <- TRUE
# }
# 
# if ((nrow(vlphen) <= 1) || (sum(!is.na(vlphen$Outcome)) <= 1)) {
#     stop("Invalid validation cohort size: N=", sum(!is.na(vlphen$Outcome)))
# }
# 
# if (binary.phen) {
#     if (sum(vlphen$Outcome == 1, na.rm=TRUE) <= 1) {
#         stop("Too few cases: N_case=", sum(vlphen$Outcome == 1, na.rm=TRUE))
#     }
# 
#     if (sum(vlphen$Outcome == 0, na.rm=TRUE) <= 1) {
#         stop("Too few controls: N_control=", sum(vlphen$Outcome == 0, na.rm=TRUE))
#     }
# }
valphenofile <- paste(tempprefix, "/", valtag, ".phen", sep='')
cat("Producing a combined prediction model...\n")

# Combined average 
# vlY <- vlphen$Outcome

# prisk <- rep(0, length(vlY))    
prisk <- rep(0, samplesize)  

for (WINSHIFT in WINSHIFT.list) {
    
    # prisk0 <- rep(0, length(vlY))
    prisk0 <- rep(0, samplesize)
    
    for (chr in 1:22) {
        ## Read per-chrom genetic risk file

        prisk.file <-
            paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                  ".predY_pg.", valtag, ".chrom", chr, ".qctoolout", sep='')

        if (file.exists(prisk.file)) {
            if (chr == 1) {
                cat("Reading", prisk.file, "(bgen)...\n")
            } else {
                cat("...", chr, "...")
            }
            
            prisk.tab <- read.table(prisk.file, header=TRUE, comment="#",
                                    stringsAsFactors=FALSE)

            # FIXME
            ASSERT(all(prisk.tab$sample == vlfam[, 1]) ||
                   all(prisk.tab$sample == vlfam[, 2]))

            prisk.chr <- prisk.tab$NPS_risk_score
            
        } else {
            prisk.file <-
                paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                      ".predY_pg.", valtag, ".chrom", chr, ".sscore", sep='')

            prisk.tab <- read.delim(prisk.file, header=FALSE, sep="\t")
            # prisk.tab <- try(read.delim(prisk.file, header=FALSE, sep="\t"), 
            #                  silent = T)
            # if (inherits(prisk.tab, "try-error") == T) {
            #     cat(paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
            #               ".predY_pg.", valtag, ".chrom", chr, ".sscore is wrong\n", sep=''))
            #     prisk.tab <- rep(0, samplesize) 
            # }
            
            if (ncol(prisk.tab) == 5) {
                ## PLINK2 generated

                if (chr == 1) {
                    cat("Reading", prisk.file, "(plink2)...\n")
                } else {
                    cat("...", chr, "...")
                }
                
                prisk.tab <- read.delim(prisk.file, header=TRUE, sep="\t")

                ASSERT(all(c("IID", "NMISS_ALLELE_CT", "SCORE1_AVG") %in%
                           colnames(prisk.tab)))
                ASSERT(colnames(prisk.tab)[2] == "IID")

                ## Reorder IID
                rownames(prisk.tab) <-
                    paste(prisk.tab[, 1], prisk.tab[, 2], sep=":")
                
                prisk.tab <-
                    prisk.tab[paste(vlfam[, 1], vlfam[, 2], sep=":"), ]
            
                ASSERT(all(prisk.tab[, 1] == vlfam[, 1]))
                ASSERT(all(prisk.tab[, 2] == vlfam[, 2]))
                ASSERT(all(!is.na(prisk.tab$SCORE1_AVG)))
                
                ## get scores
                prisk.chr <- prisk.tab$SCORE1_AVG * prisk.tab$NMISS_ALLELE_CT
                
            } else {
                if (chr == 1) {
                    cat("Reading", prisk.file, "(nps_generic)...\n")
                } else {
                    cat("...", chr, "...")
                }
                
                prisk.chr <- prisk.tab[, 1]
            }
        }

        adjbetahat.file <-
            paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                  ".adjbetahat_tail.chrom", chr, ".txt", sep='')
        adjbetahats <- read.delim(adjbetahat.file, header=FALSE)[, 1]
   
        if (any(adjbetahats != 0)) {
        
            prisk.file <-
                paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                      ".predY_tail.", valtag, ".chrom", chr, ".qctoolout",
                      sep='')
        
            if (file.exists(prisk.file)) {

                if (chr == 1) {
                    cat("Reading", prisk.file, "(bgen)...\n")
                } else {
                    cat(chr, "")
                }

                if (file.info(prisk.file)$size > 0) {
            
                    prisk.tab <-
                        read.table(prisk.file, header=TRUE, comment="#",
                                   stringsAsFactors=FALSE)

                    ASSERT(all(prisk.tab$sample == vlfam[, 1]) ||
                           all(prisk.tab$sample == vlfam[, 2]))
                
                    prisk.chr <- prisk.chr + prisk.tab$NPS_risk_score
                }
            
            } else {
                prisk.file <-
                    paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                          ".predY_tail.", valtag, ".chrom", chr, ".sscore",
                          sep='')
            
                prisk.tab <- read.delim(prisk.file, header=FALSE, sep="\t")

                prisk.tab <- try(read.delim(prisk.file, header=FALSE), silent = T)
                if(inherits(prisk.tab, "try-error") == T){
                    cat(paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                              ".predY_tail.", valtag, ".chrom", chr, ".sscore is wrong\n",
                              sep=''))
                    prisk.tab <- rep(0, samplesize)
                }

                
                if (ncol(prisk.tab) == 5) {
                    ## PLINK2 generated

                    if (chr == 1) {
                        cat("Reading", prisk.file, "(plink2)...\n")
                    } else {
                        cat("...", chr, "...")
                    }
                
                    prisk.tab <- read.delim(prisk.file, header=TRUE, sep="\t")

                    ASSERT(all(c("IID", "NMISS_ALLELE_CT", "SCORE1_AVG") %in%
                               colnames(prisk.tab)))
                    ASSERT(colnames(prisk.tab)[2] == "IID")
        
                    ## Reorder IID
                    rownames(prisk.tab) <-
                        paste(prisk.tab[, 1], prisk.tab[, 2], sep=":")
                
                    prisk.tab <-
                        prisk.tab[paste(vlfam[, 1], vlfam[, 2], sep=":"), ]
            
                    ASSERT(all(prisk.tab[, 1] == vlfam[, 1]))
                    ASSERT(all(prisk.tab[, 2] == vlfam[, 2]))
                    ASSERT(all(!is.na(prisk.tab$SCORE1_AVG)))
                
                    ## get scores
                    prisk.chr <- prisk.chr +
                        prisk.tab$SCORE1_AVG * prisk.tab$NMISS_ALLELE_CT
                
                } else {
                    if (chr == 1) {
                        cat("Reading", prisk.file, "(nps_generic)...\n")
                    } else {
                        cat("...", chr, "...")
                    }
                    
                    prisk.chr <-
                        prisk.chr + prisk.tab[, 1]
                }
            }
        }
        
        prisk0 <- prisk0 + prisk.chr
        
    }
    
    prisk <- prisk + prisk0
    
}

cat("\n")

# sex.baseline <- 0
# 
# for (WINSHIFT in WINSHIFT.list) {
#     covariate.file <-
#         paste(tempprefix, "win_", WINSHIFT, ".covariate.RDS", sep='')
#     if (file.exists(covariate.file)) {
#         sex.baseline <- sex.baseline + readRDS(covariate.file)$sex
#     } 
# }
# 
# if (sex.baseline != 0) {
#     cat("Adding sex baseline:", sex.baseline, "...\n")
# 
#     if (any((vlfam[, 5] != 1) & (vlfam[, 5] != 2))) {
#         stop("Unkown sex:", valfamfile)
#     }
#     
#     sex.covariate <- (vlfam[, 5] - 1)
# 
#     prisk <- prisk + (sex.baseline * sex.covariate)
# }


filename <- paste(valphenofile, ".nps_score", sep='')
# df.out <- cbind(vlphen, Score=prisk)
# write.table(df.out, file=filename,
#             row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(prisk, file=filename,
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
cat("Done\nPolygenic scores are saved in", filename, ".\n")
# 
# cat("Observed-scale R2 =", cor(vlY, prisk, use="pairwise")**2, "\n")
# 
# if ("TotalLiability" %in% colnames(vlphen)) {
#     ## For simulated phenotypes
#     vlL <- vlphen$TotalLiability
#     cat("Liability-scale R2 =", cor(vlL, prisk, use="pairwise")**2, "\n")
# }

# if (binary.phen) {
#     if (require(pROC, quietly=TRUE)) {
#         
#         library(pROC)
# 
#         cat("AUC:\n")
#         print(roc(cases=prisk[which(vlY == 1)], controls=prisk[which(vlY == 0)], ci=TRUE))
#         
#     } else {
#         cat("Skip AUC calculation\n")
#         cat("Please install pROC package to enable this\n")
#     }
#     
#     if (require(DescTools, quietly=TRUE)) {
#         
#         library(DescTools)
#     
#         mod <- glm(vlY ~ prisk, family=binomial(link="logit"))
#     
#         cat("Nagelkerke's R2 =", PseudoR2(mod, "Nagelkerke"), "\n")
#         
#         ## print(mod)
#         
#     } else {
#         cat("Skip Nagelkerke's R2 calculation\n")
#         cat("Please install DescTools package to enable this\n")
#     }
# 
#     ## Tail OR
#     prisk <- prisk[!is.na(vlY)]
#     vlY <- vlY[!is.na(vlY)]
# 
#     cutoff <- 0.05
#     cutoff.at <- round(length(prisk) * as.numeric(cutoff))
#     cutoff.prisk <- sort(prisk, decreasing=TRUE)[cutoff.at]
#     
#     odds1 <- sum(prisk >= cutoff.prisk & vlY == 1) /
#         sum(prisk >= cutoff.prisk & vlY == 0)
#     odds0 <- sum(prisk < cutoff.prisk & vlY == 1) /
#         sum(prisk < cutoff.prisk & vlY == 0)
# 
#     cat("Tail OR (5%):", (odds1 / odds0), "\n")
# 
#     ## calculate CI on Tail OR (time-consuming computation)
#     
#     ## library(boot)
# 
#     ## bt.iter <- 1000
# 
#     ## btOR <- function(data, ind, cutoff) {
#     ##     d <- data[ind, ]
# 
#     ##     odds1 <- sum(d$prisk >= cutoff & d$vlY == 1) /
#     ##         sum(d$prisk >= cutoff & d$vlY == 0)
#     ##     odds0 <- sum(d$prisk < cutoff & d$vlY == 1) /
#     ##         sum(d$prisk < cutoff & d$vlY == 0)
# 
#     ##     odds1 / odds0
#     ## }
# 
#     ## btdist <-
#     ##     boot(data=data.frame(prisk=prisk, vlY=vlY), statistic=btOR, R=bt.iter,
#     ##          cutoff=cutoff.prisk)
# 
#     ## print(btdist$t0)
#     ## print(boot.ci(btdist, type="perc"))
# }

cat("Done\n")

