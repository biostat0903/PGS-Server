VERSION <- "1.1"

cat("Non-Parametric Shrinkage", VERSION, "\n")

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

#########################################################################

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) < 1) {    
  stop("Usage: Rscript nps_reweight.R <work dir> [<WINSHIFT> ...]")
}

tempprefix <- paste(cargs[1], "/", sep='')
# ys
p_num <- as.numeric(cargs[2])
# ys

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))
summstatfile <- args[["summstatfile"]] 
traindir <- args[["traindir"]]
trainfreqfile <- args[["trainfreqfile"]]
traintag <- args[["traintag"]]
trainfamfile <- args[["trainfamfile"]]
trphen <- args[["trainpheno"]]
WINSZ <- args[["WINSZ"]]

ASSERT(!is.null(trphen))

# if (length(cargs) > 1) {
# 
#   WINSHIFT.list <- as.numeric(cargs[2:length(cargs)])
# 
# } else {

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
# }

if (any(is.na(WINSHIFT.list)) || any(WINSHIFT.list < 0) ||
    any(WINSHIFT.list >= WINSZ)) {

  if (length(cargs) > 1) {
    stop("Invalid shift (window size =", WINSZ, "):",
         cargs[2:length(cargs)])
  } else {
    stop("Invalid shift (window size =", WINSZ, "):",
         WINSHIFT.list)
  }
}

#########################################################################

predY <- NULL

for (WINSHIFT in WINSHIFT.list) {

  cat("----- Shifted by", WINSHIFT, "-----\n")


  ## Load partition data
  part <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".part.RDS", sep=''))

  Nt <- part[["Nt"]]
  nLambdaPT <- part[["nLambdaPT"]]
  nEtaPT <- part[["nEtaPT"]]
  lambda.q <- part[["lambda.q"]]
  betahatH.q <- part[["betahatH.q"]]

  #########################################################################

  ## phenotypes

  # ys
  trfam <- read.delim(trainfamfile, sep=" ", header=FALSE,
                      stringsAsFactors=FALSE)[, c(1, 2, 3, 4, 5, 5+p_num)]
  # ys

  if (ncol(trfam) != 6) {
    ## re-try with tab delimination

    # ys
    trfam <- read.delim(trainfamfile, sep="\t", header=FALSE,
                        stringsAsFactors=FALSE)[, c(1, 2, 3, 4, 5, 5+p_num)]
    # ys
  }

  ASSERT(ncol(trfam) == 6)
  ## ASSERT(all(!is.na(trphen$Outcome)))
  ASSERT(all(trphen$FID == trfam[, 1]))
  ASSERT(all(trphen$IID == trfam[, 2]))

  use.lda <- TRUE
  binary.phen <- TRUE

  if (length(unique(trphen$Outcome[!is.na(trphen$Outcome)])) > 4) {
    ## Quantitative phenotypes
    cat("Quantitative phenotype detected:\n")
    cat("train fam file:", trainfamfile, "\n")
    cat("Size of training cohort:", Nt, "\n")
    cat("# of missing phenotypes:", sum(is.na(trphen$Outcome)), "\n")
    cat("# of non-missing phenotypes:", sum(!is.na(trphen$Outcome)), "\n")

    binary.phen <- FALSE
    use.lda <- FALSE

  } else {
    ## Binary phenotypes
    cat("Binary phenotype detected:\n")
    cat("train fam file:", trainfamfile, "\n")
    cat("Size of training cohort:", Nt, "\n")
    cat("# of missing phenotypes:", sum(is.na(trphen$Outcome)), "\n")
    cat("# of cases:", sum(trphen$Outcome == 1, na.rm=TRUE), "\n")
    cat("# of controls:", sum(trphen$Outcome == 0, na.rm=TRUE), "\n")

    binary.phen <- TRUE
    use.lda <- TRUE
  }

  trY <- trphen$Outcome

  ASSERT(Nt == length(trY))

  ## Check covariates
  sex.covariate <- NULL

  if (any(trfam[, 5] == 1) && any(trfam[, 5] == 2)) {

    cat("Sex covariate detected.\n")

    ## Unknown sex not allowed
    if (any((trfam[, 5] != 1) & (trfam[, 5] != 2))) {
      stop("Sex is either 1 or 2; unknown sex is not allowed")
    }

    use.lda <- FALSE
    sex.covariate <- (trfam[, 5] - 1)
  }

  if (use.lda) {
    cat("Using linear discriminary analysis...\n")
  } else {
    cat("Using linear regression...\n")
  }

  #########################################################################
  ## Read partitions

  trPT <- array(0, dim=c(Nt, nLambdaPT, nEtaPT, 1))

  for (chrom in 1:22) {

    trPT.chr <-
      readRDS(paste(tempprefix, "win_", WINSHIFT, ".trPT.", chrom,
                    ".RDS", sep=''))

    trPT <- trPT + trPT.chr
  }
  
  #########################################################################

  PTwt <- array(0, dim=c(nLambdaPT, nEtaPT, 1))

  for (I in 1:nLambdaPT) {

    for (J in 1:nEtaPT) {
      K <- 1

      if (use.lda) {

        trcaVAR <- var(trPT[which(trY == 1), I, J, K])
        trctVAR <- var(trPT[which(trY == 0), I, J, K])
        trptVAR <- (trcaVAR + trctVAR) / 2

        trcaMU <- mean(trPT[which(trY == 1), I, J, K])
        trctMU <- mean(trPT[which(trY == 0), I, J, K])

        PTwt[I, J, K] <- (trcaMU - trctMU) / trptVAR

      } else {
        ## Use linear regression
        x <- trPT[, I, J, K]

        if (is.null(sex.covariate)) {
          trlm <- lm(trY ~ x)
        } else {
          trlm <- lm(trY ~ x + sex.covariate)
        }

        PTwt[I, J, K] <- trlm$coefficients[2]
      }

      if (is.nan(PTwt[I, J, K])) {
        cat("WARNING: partition (", I, ",", J, ") produced NaN\n")
        PTwt[I, J, K] <- 0
      }

    }
  }

  predY0 <- rep(0, Nt)

  for (I in 1:nLambdaPT) {

    for (J in 1:nEtaPT) {
      K <- 1

      predY0 <- predY0 + PTwt[I, J, K] * trPT[, I, J, K]
    }
  }
  cat(summary(trPT[, I, J, K]), "\n")
  # cat(PTwt[ , , 1])

  cat("Saving ", nLambdaPT, "x", nEtaPT, "partition weights...")

  saveRDS(PTwt, paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS", sep=''))

  cat("OK\n")

  #########################################################################

  ## tail partition
  trPT.tail <- rep(0, Nt)

  for (chrom in 1:22) {

    trPT.tail.file <-
      paste(tempprefix, "trPT.", chrom, ".tail.RDS", sep='')

    # FIXME
    #        ASSERT(file.exists(trPT.tail.file))

    if (file.exists(trPT.tail.file)) {
      cat("Loading S0 partition for chrom", chrom, "...\n")

      trPT.tail.chr <- readRDS(trPT.tail.file)

    } else {
      cat("No GWAS-sig S0 partition for chrom", chrom, ": set to 0...\n")

      trPT.tail.chr <- rep(0, Nt)
    }

    trPT.tail <- trPT.tail + trPT.tail.chr
  }
  
  PTwt.tail <- 0
  sex.baseline <- 0
  
  if (any(trPT.tail != 0)) {

    if (binary.phen) {
      
      ## Use logistic regression
      if (is.null(sex.covariate)) {

        trlm <-
          glm(trY ~ predY0 + trPT.tail,
              family=binomial(link="logit"))

      } else {

        trlm <-
          glm(trY ~ predY0 + trPT.tail + sex.covariate,
              family=binomial(link="logit"))

        sex.baseline <- trlm$coefficients[4] / trlm$coefficients[2]
      }

    } else {
      
      ## Use linear regression
      if (is.null(sex.covariate)) {

        trlm <- lm(trY ~ predY0 + trPT.tail)

      } else {

        trlm <- lm(trY ~ predY0 + trPT.tail + sex.covariate)

        sex.baseline <- trlm$coefficients[4] / trlm$coefficients[2]
      }
    }

    PTwt.tail <- trlm$coefficients[3] / trlm$coefficients[2]

    ## tail.h2.ratio <-
    ##    (cor(trY, as.vector(trPT.tail * trlm$coefficients[3]), use="pairwise")**2) /
    ##    (cor(trY, as.vector(predY0 * trlm$coefficients[2]), use="pairwise")**2)
    ## print(tail.h2.ratio)

  } else {
    ## No S0 partition
    if (binary.phen) {
      ## Use logistic regression
      trlm <-
        glm(trY ~ predY0 + sex.covariate,
            family=binomial(link="logit"))
    } else {
      ## Use linear regression
      cat(summary(predY0), "\n")
      trlm <- lm(trY ~ predY0 + sex.covariate)
    }

    sex.baseline <- trlm$coefficients[3] / trlm$coefficients[2]

  }

  cat("Weight for S0 =", PTwt.tail, "\n")

  cat("Saving S0 weight...")

  saveRDS(PTwt.tail,
          paste(tempprefix, "win_", WINSHIFT, ".PTwt.tail.RDS", sep=''))

  cat("OK\n")

  if (!is.null(sex.covariate)) {
    cat("Female baseline =", sex.baseline, "\n")

    saveRDS(list(sex=sex.baseline),
            paste(tempprefix, "win_", WINSHIFT, ".covariate.RDS", sep=''))
  }

  ######################################################################
  ## Training R2

  #    cat("Observed scale R2 in training =",
  #        cor(trY, predY0 + PTwt.tail * trPT.tail)**2, "\n")

  if (is.null(predY)) {
    predY <- predY0 + PTwt.tail * trPT.tail
  } else {
    predY <- predY + (predY0 + PTwt.tail * trPT.tail)
  }

  if (!is.null(sex.covariate)) {
    predY <- predY + (sex.baseline * sex.covariate)
  }

  #########################################################################
  ## back2snpeff

  for (CHR in 1:22) {

    cat("Re-weighting SNPs in chr", CHR, "...\n")

    ## Read summary stats (discovery)
    summstat.chr <- read.delim(paste(summstatfile, ".", CHR, sep=''),
                               header=TRUE, stringsAsFactors=FALSE,
                               sep="\t")
    #dim(summstat)

    ## Use traing AF instead of discovery AF
    trfrq.chr <-
      read.table(paste(trainfreqfile, ".", CHR, sep=''), header=TRUE)
    tr.se.chr <- sqrt(2 * trfrq.chr$AAF * (1 - trfrq.chr$AAF))
    #plot(tr.se, se, cex=0.25)
    #abline(0, 1, col="red")

    M.chr <- length(tr.se.chr)

    ASSERT(M.chr == nrow(summstat.chr))

    cat("M", "CHR", CHR, "=", M.chr, "\n")

    wt.betahat <- c()

    I <- 1

    winfilepre <-
      paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')

    while (file.exists(paste(winfilepre, ".pruned", ".table", sep=''))) {

      tailfixfile <- paste(winfilepre, ".pruned", ".table", sep='')

      wintab <- read.delim(tailfixfile, header=TRUE, sep="\t")

      lambda0 <- wintab$lambda
      etahat0 <- wintab$etahat

      Nq <- sum(lambda0 > lambda.q[1])

      if (Nq == 0) {
        ## No projection left

        ## Fill in with 0
        windata <- readRDS(file=paste(winfilepre, ".RDS", sep=''))
        s0 <- windata[["eigen"]]
        M <- nrow(s0$vectors)   # number of SNPs in the window
        wt.betahat <- c(wt.betahat, rep(0, M))

        ## move on to next iteration
        I <- I + 1

        winfilepre <-
          paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I,
                sep='')

        next
      }

      Q0 <- readRDS(paste(winfilepre, ".Q.RDS", sep=''))

      etahat0 <- etahat0[lambda0 > lambda.q[1]]
      Q0 <- Q0[, lambda0 > lambda.q[1], drop=FALSE]
      lambda0 <- lambda0[lambda0 > lambda.q[1]]

      wt0 <- rep(NA, Nq)

      for (Il in 1:nLambdaPT) {

        lambda.lo <- lambda.q[Il]
        lambda.hi <- lambda.q[Il+1]
        in.lambda.bin <- lambda0 > lambda.lo & lambda0 <= lambda.hi

        for (Je in 1:nEtaPT) {

          betahatH.lo <- betahatH.q[Je, Il]
          betahatH.hi <- betahatH.q[Je+1, Il]
          in.betahatH.bin <-
            (in.lambda.bin &
               abs(etahat0) > betahatH.lo &
               abs(etahat0) <= betahatH.hi)

          if (any(in.betahatH.bin)) {
            wt0[in.betahatH.bin] <- PTwt[Il, Je, 1]
          }
        }
      }

      if (any(etahat0 == 0)) {
        wt0[etahat0 == 0] <- 0
      }

      ASSERT(all(!is.na(wt0)))

      #       Compared to manuscript, we did not scale qX0 with lambda^(-1/2),
      #       thus no need to scale here again, wt0 includes the factor already.
      #       etahat0.adj <- etahat0 * wt0 / sqrt(lambda0)
      #       In effect, however, this is equivalent to the description of methods
      #       as this is compensated at the later step.
      #       The difference is intended to minimize numerical errors.
      etahat0.adj <- etahat0 * wt0

      wt.betahat <- c(wt.betahat, Q0 %*% as.matrix(etahat0.adj))

      ASSERT(all(!is.na(wt.betahat)))

      ## move on to next iteration
      I <- I + 1

      winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')
    }

    ## pad
    M.written <- length(wt.betahat)

    if ((M.chr - M.written) > 0) {
      cat("Pad ", (M.chr - M.written),
          " SNPs with 0 at the end of chrom\n")

      wt.betahat <- c(wt.betahat, rep(0, M.chr - M.written))
    }

    ASSERT(M.chr == length(wt.betahat))


    ## Add tail betahats
    tailbetahatfile <- paste(tempprefix, "tail_betahat.", CHR, ".table",
                             sep='')

    ASSERT(file.exists(tailbetahatfile))

    betahat.tail.chr <-
      read.delim(tailbetahatfile, header=FALSE, sep="\t")[, 1]

    ASSERT(length(betahat.tail.chr) == M.chr)

    wt.betahat.tail <- betahat.tail.chr * PTwt.tail

    ## Account for MAF
    ## se: discovery af
    ##    wt.betahat <- wt.betahat / se[snpIdx0 + c(1:M.chr)]

    ASSERT(length(tr.se.chr) == M.chr)

    ## se: training af
    wt.betahat <- wt.betahat / tr.se.chr
    wt.betahat.tail <- wt.betahat.tail / tr.se.chr

    ## Save
    filename.pg <- paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                         ".adjbetahat_pg.chrom", CHR, ".txt", sep='')
    ## cat("Saving reweighted snpeffs:", filename.pg, "...")
    write.table(data.frame(betahat=wt.betahat),
                file=filename.pg,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    ## cat("OK\n")

    filename.tail <- paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                           ".adjbetahat_tail.chrom", CHR, ".txt", sep='')
    ## cat("Saving reweighted snpeffs:", filename.tail, "...")
    write.table(data.frame(betahat=wt.betahat.tail),
                file=filename.tail,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    ## cat("OK\n")

  }
}

cat("Observed scale R2 in training cohort =", cor(trY, predY, use="pairwise")**2, "\n")

if (require(pROC, quietly=TRUE) &&
    (length(unique(trphen$Outcome[!is.na(trphen$Outcome)])) <= 4)) {
  ## Binary phenotypes

  library(pROC)

  cat("AUC in training cohort:\n")
  print(roc(cases=predY[which(trY == 1)], controls=predY[which(trY == 0)], ci=TRUE))

}

cat("Done\n")
