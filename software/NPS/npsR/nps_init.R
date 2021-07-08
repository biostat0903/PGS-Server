VERSION <- "1.1"

cat("Non-Parametric Shrinkage", VERSION, "\n")

print.usage <- function() {
    cat("Usage:\n")
    cat("Rscript nps_init.R <summary stats file> <train dir> <train fam file> <train phenotype file> <train dataset ID> <window size> <work dir>\n")
    cat("\nRscript nps_init.R [options]\n")
    cat("    --gwas <summary stats file>\n")
    cat("    --train-dir <train dir>\n")
    cat("    --train-dataset <train dataset ID>\n")            
    cat("    --out <work dir>\n")
    cat("    --window-size <window size>\n")
    cat("          default value: 4000\n")
    cat("    --train-fam <train fam file>\n")
    cat("          default value: <train dir>/<train dataset ID>.fam\n")
    cat("    --train-phen <train phenotype file>\n")
    cat("          default value: Use <train fam file> or\n")
    cat("          <train dir>/<train dataset ID>.phen\n")
    ##ys
    cat("    --p <train phenotype number>\n")
    cat("          default value: 1\n")
    ##ys
    cat("\n\n")
}

ASSERT <- function(test) {
    if (length(test) == 0) {
        stop(paste("ASSERT fail for empty conditional:",
                   deparse(substitute(test))))
    }

    if (is.na(test)) {
        stop(paste("ASSERT fail for missing value:",
                   deparse(substitute(test))))
    }
    
    if (!test) {
        stop(paste("ASSERT fail:", deparse(substitute(test))))
    }
}


args <- commandArgs(trailingOnly=TRUE)

if (any(startsWith(args, "--"))) {
    ## Parse long options

    summstatfile <- NULL
    traindir <- NULL
    trainfamfile <- NULL
    trainphenofile <- NULL
    traintag <- NULL
    WINSZ <- NULL
    tempprefix <- NULL
    ##ys
    p <- NULL
    ##ys
    arg.index <- 1

    while (arg.index <= length(args)) {

        if (!startsWith(args[arg.index], "--")) {
            cat("Cannot recognize the option:", args[arg.index], "\n")
            print.usage()
            q()
        }

        opt <- substr(args[arg.index], 3, nchar(args[arg.index]))

        if (opt == "gwas") {
            arg.index <- arg.index + 1

            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            summstatfile <- args[arg.index]

        } else if (opt == "train-dir") {
            arg.index <- arg.index + 1

            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            traindir <- args[arg.index]

        } else if (opt == "train-dataset") {
            arg.index <- arg.index + 1

            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            traintag <- args[arg.index]
            
        } else if (opt == "out") {
            arg.index <- arg.index + 1

            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            tempprefix <- paste(args[arg.index], "/", sep='')

        } else if (opt == "window-size") {
            arg.index <- arg.index + 1

            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            WINSZ <- as.numeric(args[arg.index])

            if (is.na(WINSZ)) {
                stop("Invalid window size: ", args[arg.index])
            }

        } else if (opt == "train-fam") {
            arg.index <- arg.index + 1
            
            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            trainfamfile <- args[arg.index]
            
        } else if (opt == "train-phen") {
            arg.index <- arg.index + 1

            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }

            trainphenofile <- args[arg.index]      

        } 
        ## ys
        else if (opt == "p") {
            arg.index <- arg.index + 1
            
            if (arg.index > length(args)) {
                stop("missing option value for ", paste("--", opt, sep=''))
            }
            
            p <- as.numeric(args[arg.index])      
            
        } 
        ##ys
        else {
            cat("Cannot recognize the option: ", args[arg.index], "\n")
            print.usage()
            q()
        }
        
        arg.index <- arg.index + 1
    }

    ## Check required values
    if (is.null(summstatfile)) {
        stop("--gwas required")
    }

    if (is.null(traindir)) {
        stop("--train-dir required")
    }

    if (is.null(traintag)) {
        stop("--train-dataset required")
    }

    if (is.null(tempprefix)) {
        stop("--out required")
    }

    ## Set default values
    if (is.null(WINSZ)) {
        WINSZ <- 4000
        cat("Using the default window size:", WINSZ, "\n")
    }

    if (is.null(trainfamfile)) {
        trainfamfile <- paste(traindir, "/", traintag, ".fam", sep='')
        cat("Using the default training fam file path:", trainfamfile, "\n")
    }
    
} else {
    ## Sequential

    if (length(args) != 7) {
        print.usage()
        q()
    }
    
    summstatfile <- args[1]
    traindir <- args[2]
    trainfamfile <- args[3]
    trainphenofile <- args[4]
    traintag <- args[5]
    WINSZ <- as.numeric(args[6])
    ##ys
    p <- as.numeric(args[7])
    ##ys

    if (is.na(WINSZ)) {
        stop("Invalid window size:", args[6])
    }
    
    tempprefix <- paste(args[7], "/", sep='')    
}

if (WINSZ <= 10) {
    stop("Too small window size: ", WINSZ)
}

#################################################################
# SANITY CHECKS

if (!file.exists(summstatfile)) {
    stop("File does not exists:", summstatfile)
}

if (!dir.exists(traindir)) {
   stop("Directory does not exists:", traindir)
}

if (!file.exists(trainfamfile)) {
   stop("File does not exists:", trainfamfile)
}

summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
# dim(summstat)

if (length(intersect(colnames(summstat), 
                     c("chr", "pos", "ref", "alt", "reffreq", "pval",
                       "effalt"))) 
    != 7) {
    cat(paste(colnames(summstat), collapse="\t"), "\n")
    cat("Expected essential columns: ",
        paste(c("chr", "pos", "ref", "alt", "reffreq", "pval", "effalt"),
              collapse="\t"), "\n")
    stop("Missing essential columns:", summstatfile)
}

if (nchar(summstat$chr[1]) < 4 || substr(summstat$chr[1], 1, 3) != "chr") {
    stop("Chromosome names are expected to be chr1, ..., chr22:", summstatfile)
}

## Sort
summstat.chr <- substr(summstat$chr, 4, nchar(summstat$chr))

if (length(setdiff(summstat.chr, as.character(1:22))) != 0) {
    cat(paste(setdiff(summstat.chr, as.character(1:22)), collapse=", "), "\n")
    stop("expect only chr1 through chr22:", summstatfile)
}

summstat.chr <- as.numeric(summstat.chr)
summstat <- cbind(summstat, chrnum=summstat.chr)

summstat <- summstat[order(summstat$chrnum, summstat$pos), ]
summstat.chr <- summstat$chrnum

## if (is.unsorted(summstat.chr)) {
##     stop("Not sorted by chr:", summstatfile)
## }

## for (CHR in 1:22) {
##     if (is.unsorted(summstat$pos[summstat.chr == CHR])) {
##         stop("chr", CHR, ": not sorted by pos:", summstatfile)
##     }

##     if (any(duplicated(summstat$pos[summstat.chr == CHR]))) {
##         summstat0 <- summstat[summstat.chr == CHR, ]
##         dup.pos <- duplicated(summstat0$pos)
##         print(head(summstat0[summstat0$pos %in% dup.pos, ]))
            
##         stop("chr", CHR, ": duplicated pos or tri-allelic:", summstatfile)
##     }
## }

if (any(is.na(summstat$reffreq))) {
    stop("NA in reffreq:", summstatfile)
}

if (any(is.na(summstat$effalt))) {
    stop("NA in effalt:", summstatfile)
}

if (any(is.na(summstat$pval))) {
    stop("NA in pval:", summstatfile)
}

if (any(summstat$pval == 0)) {
    stop("pval underflow (pval == 0):", summstatfile)
}

## if (any(nchar(summstat$ref) != 1) || any(nchar(summstat$alt) != 1)) {
##     stop("InDels are not allowed:", summstatfile)
## }

summstat.SNPID <- paste(summstat.chr, ":", summstat$pos,
                        "_", summstat$ref, "_", summstat$alt,
                        sep='')

if (any(duplicated(summstat.SNPID))) {
    cat(paste(summstat.SNPID[duplicated(summstat.SNPID)], collapse=", "), "\n")
    stop("Duplicated SNPs in GWAS summary statistics file")
}

if (!dir.exists(tempprefix)) {
    dir.create(tempprefix)
}

if (!dir.exists(paste(tempprefix, "/log", sep=''))) {
    dir.create(paste(tempprefix, "/log", sep=''))
} else {
    log.files <- paste(tempprefix, "/log", "/nps_*.Rout.*", sep='')
    cat("Removing log files: ", log.files, "...")
    unlink(log.files)
    cat(" OK\n")
}


# train allele freq
trfrq.combined <- NULL

trainfreqfile <-
    paste(tempprefix, "/", traintag, ".meandos", sep='')

# snp info
trSNPID <- c()

for (CHR in 1:22) {

    # Read train allele freq
    meandosfile <-
        paste(traindir, "/chrom", CHR, ".", traintag, ".meandos", sep='')

    trfrq <- read.table(meandosfile, header=TRUE, stringsAsFactors=FALSE)
    
    if (length(colnames(trfrq)) != 2 ||
        any(colnames(trfrq) != c("SNPID", "AAF"))) {

        stop("Invalid .meandos header: ", meandosfile)
    }

    # Save AAF
    trainfreqfile.chr <- paste(trainfreqfile, ".", CHR, sep='')
        
    ## cat("Copying training AAF file to: ", trainfreqfile.chr, "...")

    write.table(trfrq, file=trainfreqfile.chr,
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    ## cat(" OK\n")

    trfrq.combined <- rbind(trfrq.combined, trfrq)

    # Read SNP infos
    snpinfofile <-
        paste(traindir, "/chrom", CHR, ".", traintag, ".snpinfo", sep='')

    snpinfo <- read.table(snpinfofile, header=TRUE, stringsAsFactors=FALSE)

    trSNPID.chr <- paste(snpinfo$chromosome, ":", snpinfo$position,
                         "_", snpinfo$alleleA, "_", snpinfo$alleleB,
                         sep='')

    if (any(duplicated(trSNPID.chr))) {
        cat(paste(trSNPID.chr[duplicated(trSNPID.chr)], collapse=", "), "\n")
        stop("Duplicated SNPs in training genotype files")
    }

    if (any(duplicated(snpinfo$position))) {
        cat("chr", CHR, ":\n")
        cat(paste(snpinfo$position[duplicated(snpinfo$position)],
                  collapse=", "), "\n")
        stop("Only biallelic SNPs are allowed in training genotypes")
             
    }
    
    if (any(!(trSNPID.chr %in% summstat.SNPID))) {
        cat(paste(setdiff(trSNPID.chr, summstat.SNPID), collapse="\n"), "\n")
        stop("Missing summary statistics for the above SNPs in training cohort")
    }

    trSNPID <- c(trSNPID, trSNPID.chr)
}

summstat <- summstat[summstat.SNPID %in% trSNPID, ]
summstat.SNPID <- summstat.SNPID[summstat.SNPID %in% trSNPID]

if (nrow(trfrq.combined) != nrow(summstat)) {
    print(nrow(summstat))
    print(nrow(trfrq.combined))
    stop("The number of markers does not match:", summstatfile, ", ", 
         paste(traindir, "/chromXX.", traintag, ".meandos", sep=''))
}

if (length(trSNPID) != nrow(summstat)) {
    print(nrow(summstat))
    print(length(trSNPID))
    stop("The number of markers does not match:", summstatfile, ", ", 
         paste(traindir, "/chromXX.", traintag, ".snpinfo", sep=''))
}

if (any(trSNPID != summstat.SNPID)) {
    cat(head(data.frame(summstat=summstat.SNPID, train=trSNPID,
                        stringsAsFactors=FALSE)[trSNPID != summstat.SNPID, ]))
    
    stop("Alleles does not align:", summstatfile, ", ",
         paste(traindir, "/chromXX.", traintag, ".snpinfo", sep=''))
}

trfam <- read.delim(trainfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)[,c(1, 2, 3, 4, 5, 5+p)]

if (ncol(trfam) != 6) {
    # re-try with tab delimination

    trfam <- read.delim(trainfamfile, sep="\t", header=FALSE,
                        stringsAsFactors=FALSE)
}

if (ncol(trfam) != 6) {    
    stop("Space or tab-delimited 6-column FAM format expected:", trainfamfile)
}

Nt <- nrow(trfam)

if (Nt <= 1) {
    stop("Invalid training cohort size:", Nt)
}

if (Nt < 1000) {
    cat("Warning: Training cohort may be too small: ", Nt, "\n")
}

trphen <- NULL

if (is.null(trainphenofile)) {
    if (all(trfam[, 6] == 0) || all(trfam[, 6] == -9)) {
        cat("Phenotype data are empty in", trainfamfile, ".\n")

        ## default phenotype file path
        trainphenofile <- paste(traindir, "/", traintag, ".phen", sep='')
        cat("Using the default training phen file path:", trainphenofile, "\n")

    } else {
        cat("Use phenotypes in", trainfamfile, ".\n")
        trphen <- trfam[, c(1, 2, 6)]
        colnames(trphen) <- c("FID", "IID", "Outcome")
    }
}

if (is.null(trphen)) {
    cat("Reading", trainphenofile, "for phenotype data.\n")

    if (!file.exists(trainphenofile)) {
        stop("File does not exists:", trainphenofile)
    }

    trphen <- read.delim(trainphenofile, sep="\t", header=TRUE,
                         stringsAsFactors=FALSE)
}

if (length(intersect(colnames(trphen), 
                     c("FID", "IID", "Outcome"))) 
    != 3) {

    stop("FID\tIID\tOutcome columns expected (tab-delimited):", trainphenfile)
}

rownames(trphen) <- paste(trphen$FID, trphen$IID, sep=":")

## missing.famIDs <- 
##     setdiff(paste(trfam[, 1], trfam[, 2], sep="\t"),
##             paste(trphen$FID, trphen$IID, sep="\t"))

## if (length(missing.famIDs) > 0) {
##     cat(paste(missing.famIDs, collapse="\n"))
##     stop("Missing phenotypes for samples in ", trainfamfile, ":",
##          trainphenofile)
## }

## if (any(is.na(trphen$Outcome))) {
##     stop("NA is not allowed for training phenotypes")
## }

if (!is.numeric(trphen$Outcome)) {
    stop("phenotype values are not numeric")
}


if (length(unique(trphen$Outcome[!is.na(trphen$Outcome)])) > 4) {
    # Quantitative phenotypes
    cat("Quantitative phenotype: Outcome\n")

    if (any(trphen$Outcome == -9)) {
        cat("WARNING: Phenotype value of \"-9\" will be interpreted as NA\n")
    }

    trphen$Outcome[which(trphen$Outcome == -9)] <- NA

} else {
    # Binary phenotypes    
    cat("Binary phenotype: Outcome\n")

    if (all(trphen$Outcome == 0 | trphen$Outcome == 1 |
            trphen$Outcome == -9, na.rm=TRUE)) {
        ## NPS v1.0 used 0/1/-9 encoding

        ## if (any(trphen$Outcome == -9)) {
        ##     stop("NA (\"-9\") is not allowed for training phenotypes")
        ## }

        trphen$Outcome[which(trphen$Outcome == -9)] <- NA
        
    } else {
        ## 1/2/0/-9 encoding
        if (!all(trphen$Outcome == 1 | trphen$Outcome == 2 |
                 trphen$Outcome == -9 | trphen$Outcome == 0, na.rm=TRUE)) {
            stop("Binary trait values should be 1, 2, 0 or -9")
        }

        if (sum(trphen$Outcome == 1, na.rm=TRUE) == 0) {
            stop("There is no control sample (\"1\") in the training data")
        }
        
        if (sum(trphen$Outcome == 2, na.rm=TRUE) == 0) {
            stop("There is no case sample (\"2\") in the training data")
        }        
        
        ## for backward compatibility
        ## Outcome code 2 -> 1 (case)
        ## Outcome code 1 -> 0 (control)
        ## Outcome code 0 -> -9 (missing)
        ## Outcome code -9 -> -9 (missing)

        ## if (any(trphen$Outcome == -9)) {
        ##     stop("NA (\"-9\") is not allowed for training phenotypes")
        ## }
        ## if (any(trphen$Outcome == 0)) {
        ##     stop("NA (\"0\") is not allowed for training phenotypes")
        ## }

        trphen$Outcome[which(trphen$Outcome == -9)] <- NA
        trphen$Outcome[which(trphen$Outcome == 0)] <- NA
        
        ## recode to 0/1
        trphen$Outcome <- trphen$Outcome - 1
    }

    if (sum(trphen$Outcome == 0, na.rm=TRUE) == 0) {
        stop("There is no control sample in the training data")
    }

    if (sum(trphen$Outcome == 1, na.rm=TRUE) == 0) {
        stop("There is no case sample in the training data")
    }
}

trphen <- trphen[paste(trfam[, 1], trfam[, 2], sep=":"), ]
trphen$FID <- trfam[, 1]
trphen$IID <- trfam[, 2]

## ASSERT(all(!is.na(trphen$Outcome)))
## ASSERT(all(trphen$FID == trfam[, 1]))
## ASSERT(all(trphen$IID == trfam[, 2]))

##################################################################
# SAVE

# Save summary statistics
summstatfile2 <- paste(tempprefix, "/harmonized.summstats.txt", sep='')

cat("Dumping harmonized summary stats file: ", summstatfile2, "...")

cols <- c("chr", "pos", "ref", "alt", "reffreq", "pval", "effalt")

if ("N" %in% colnames(summstat)) {
   cols <- c(cols, "N")
}

write.table(summstat[, cols],
            file=summstatfile2, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
cat(" OK\n")

# Save config
cat("Writing config file ...\n")

args <- list()

args[["VERSION"]] <- VERSION
args[["summstatfile"]] <- summstatfile2
args[["Nt"]] <- Nt
args[["traindir"]] <- traindir
args[["trainfamfile"]] <- trainfamfile
args[["trainfreqfile"]] <- trainfreqfile
args[["trainpheno"]] <- trainphenofile
args[["traintag"]] <- traintag
args[["WINSZ"]] <- WINSZ
args[["CXWCOR.CO"]] <- 0.3

print(args)

args[["trainpheno"]] <- trphen

saveRDS(args, file=paste(tempprefix, "args.RDS", sep=''))

cat("Done\n")
