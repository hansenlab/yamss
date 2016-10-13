## Steps
## 1. Read raw data
## 2. Background correction
## 3. Retention time alignment
## 4. Density estimation
## 5. Get cutoffs for density estimation
## 6. Get blobs
## 7. Get XICs
## 8. Quantify
## 9. Differential analysis

backgroundCorrection <- function(object, verbose = FALSE) {
    stopifnot(is(object, "CMSraw"))
    rawDT <- .rawDT(object)
    setkey(rawDT, mz, scan, sample)
    mzbreaks <- unique(c(seq(.minMZ(object), .maxMZ(object), by = 10), .maxMZ(object)))
    scanbreaks <- seq(1, .maxScan(object), 40)
    scanbreaks[length(scanbreaks)] <- .maxScan(object)
    if(verbose) {
        message("[backgroundCorrection] Get marginal intensities")
    }
    ptime1 <- proc.time()
    densGrid <- lapply(seq_len(length(mzbreaks)-1), function(m) {
        mzseq <- seq(as.integer(mzbreaks[m]*1e5), as.integer(mzbreaks[m+1]*1e5))
        DTmz <- rawDT[.(mzseq), nomatch = 0]
        lapply(seq_len(length(scanbreaks)-1), function(rt) {
            DTmzscan <- DTmz[scan %in% scanbreaks[rt]:scanbreaks[rt+1]]
            lapply(.sampleNumber(object), function(s) {
                logintens <- log2(DTmzscan[sample==s, intensity] + 1)
                if (length(logintens) < 10)
                    return(NA)
                density(logintens)
            })
        })
    })
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[backgroundCorrection] Get marginal intensities .. done in %.1f secs.", stime))
    }
    ## Estimate retention time window-specific background levels
    bgsd <- 0 # SD of normal distribution characterizing noise intensities
    r <- dnorm(bgsd)/dnorm(0)
    bgmeans <- lapply(seq_len(length(mzbreaks)-1), function(m) {
        lapply(seq_len(length(scanbreaks)-1), function(rt) {
            lapply(.sampleNumber(object), function(s) {
                dens <- densGrid[[m]][[rt]][[s]]
                if (class(dens)=="density") {
                    indexFirstPeak <- which.max(diff(dens$y) < 0)
                    densCutoff <- dens$y[indexFirstPeak]*r
                    bgindex <- which.max(dens$y <= densCutoff & c(rep(FALSE, indexFirstPeak), rep(TRUE, length(dens$y)-indexFirstPeak)))
                    return(dens$x[bgindex])
                }
                return(NA)
            })
        })
    })
    if(verbose) {
        message("[backgroundCorrection] Get region-specific background trends")
    }
    ptime1 <- proc.time()
    smooths <- lapply(.sampleNumber(object), function(s) {
        ## rows = scans, cols = M/Z bins
        ## Each col is the background trend across scans for a particular M/Z region
        bgmeanmatThisSample <- do.call(cbind, lapply(bgmeans, function(mzList) {
                                                  sapply(mzList, function(scanList) { scanList[[s]] })
                                              }))
        keepcols <- colSums(!is.na(bgmeanmatThisSample))!=0
        bgmeanmatThisSample <- bgmeanmatThisSample[,keepcols, drop = FALSE]
        ## Each col is a smoothed background trend across scans for a particular M/Z region
        bgmeanmatSmoothed <- do.call(cbind, lapply(seq_len(ncol(bgmeanmatThisSample)), function(i) {
                                                dists <- (1:ncol(bgmeanmatThisSample))-i
                                                wts <- dnorm(dists/4)
                                                weightmat <- matrix(wts, nrow = nrow(bgmeanmatThisSample),
                                                                    ncol = ncol(bgmeanmatThisSample), byrow = TRUE)
                                                df <- data.frame(intens = as.numeric(bgmeanmatThisSample),
                                                                 scan = rep(head(scanbreaks, -1),
                                                                            times = ncol(bgmeanmatThisSample)),
                                                                 weight = as.numeric(weightmat))
                                                df <- df[complete.cases(df),]
                                                lofit <- loess(intens ~ scan, data = df, weights = weight, span = 0.1)
                                                predict(lofit, .minScan(object):.maxScan(object))
                                            }))
        mzbounds <- cbind(head(mzbreaks, -1), tail(mzbreaks, -1))
        mzbounds <- mzbounds[keepcols,,drop = FALSE]
        attr(bgmeanmatSmoothed, "mzbounds") <- mzbounds
        return(bgmeanmatSmoothed)
    })
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[backgroundCorrection] Get region-specific background trends .. done in %.1f secs.", stime))
    }

    ## Perform background correction
    setkey(rawDT, mz, sample)
    if(verbose) {
        message("[backgroundCorrection] Correct intensities")
    }
    ptime1 <- proc.time()
    bgcorrDT <- rbindlist(lapply(.sampleNumber(object), function(s) {
        bgmeanmatSmoothed <- smooths[[s]]
        mzbounds <- attr(bgmeanmatSmoothed, "mzbounds")
        rbindlist(lapply(seq_len(nrow(mzbounds)), function(i) {
            mzseq <- seq(as.integer(mzbounds[i,1]*1e5), as.integer(mzbounds[i,2]*1e5))
            bgtrend <- bgmeanmatSmoothed[,i]
            bgtrend[is.na(bgtrend)] <- 0
            rawDT[.(mzseq), nomatch = 0][sample==s][,bg := (2^bgtrend[scan])-1]
        }))
    }))
    bgcorrDT[, intensity := intensity - bg]
    bgcorrDT[, bg := NULL]
    bgcorrDT <- bgcorrDT[intensity > 0]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[backgroundCorrection] Correct intensities .. done in %.1f secs.", stime))
    }
    out <- as(object, "CMSproc")
    .bgcorrDT(out) <- bgcorrDT
    return(out)
}

rtAlignment <- function(object, verbose = FALSE) {
    stopifnot(is(object, "CMSproc"))
    rawDT <- .rawDT(object)
    bgcorrDT <- .bgcorrDT(object)
    if(verbose) {
        message("[rtAlignment] Get rough M/Z regions to align")
    }
    ## Get density of M/Z values
    ## This density will be thresholded to yield M/Z regions
    mzdens <- density(bgcorrDT[,mz]/1e5, weights = bgcorrDT[,intensity]/sum(bgcorrDT[,intensity]), n = 2^ceiling(log2(nrow(bgcorrDT))), bw = 0.005)
    ## Get M/Z regions for a variety of quantile cutoffs
    qdens <- quantile(mzdens$y, seq(0.5,0.99,0.01))
    mzbounds <- lapply(qdens, function(cutoff) {
        ir <- whichAsIRanges(mzdens$y > cutoff)
        cbind(mzdens$x[start(ir)], mzdens$x[end(ir)])
    })
    ## Get 90th percentile of M/Z widths for each set of M/Z bounds
    p90 <- sapply(mzbounds, function(mat) {
        quantile(mat[,2] - mat[,1], 0.9)
    })
    ## What is the first cutoff index that for which the 90th percentile of M/Z widths
    ## is less than 0.05?
    ## i.e. We want the lowest cutoff such that the wide M/Z widths are not too wide
    wh <- which.max(p90 < 0.05)
    irmzr <- IRanges(start = as.integer(mzbounds[[wh]][,1]*1e5), end = as.integer(mzbounds[[wh]][,2]*1e5))
    if(verbose) {
        message("[rtAlignment] Get EICs for these regions")
    }
    ptime1 <- proc.time()
    eics <- getEICS(object, mzranges = irmzr)
    scans <- .minScan(object):.maxScan(object)
    eicsImputed <- lapply(eics, function(x) {
        do.call(cbind, lapply(seq_len(ncol(x)), function(col) {
                           bool <- x[,col] > 1e-6
                           if (sum(bool) < 2)
                               return(x[,col])
                           approx(scans[bool], x[bool,col], xout = scans, rule = 2)$y
                       }))
    })
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[rtAlignment] Get EICs for these regions .. done in %.1f secs.", stime))
        message("[rtAlignment] Find best shifts")
    }
    ptime1 <- proc.time()
    shifts <- -20:20
    bgcorrDT[, scanorig := scan]
    rawDT[, scanorig := scan]
    shiftsList <- lapply(seq_along(irmzr), function(i) {
        eicmat <- eics[[i]]
        eicimpmat <- eicsImputed[[i]]
        refsamp <- which.max(colSums(eicmat))
        bestShiftBySample <- sapply(.sampleNumber(object), function(s) {
            if (s==refsamp) {
                return(0)
            }
            corrShifts <- sapply(shifts, function(shift) {
                if (shift < 0) {
                    x <- tail(eicimpmat[,s], shift)
                    ref <- head(eicimpmat[,refsamp], shift)
                } else if (shift==0) {
                    x <- eicimpmat[,s]
                    ref <- eicimpmat[,refsamp]
                } else {
                    x <- head(eicimpmat[,s], -shift)
                    ref <- tail(eicimpmat[,refsamp], -shift)
                }
                cor(x, ref)
            })
            if (sum(!is.na(corrShifts))==0) {
                return(0)
            }
            shifts[which.max(corrShifts)]
        })
        return(bestShiftBySample)
    })
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[rtAlignment] Find best shifts .. done in %.1f secs.", stime))
        message("[rtAlignment] Remap scans")
    }
    ptime1 <- proc.time()
    setkey(bgcorrDT, mz, sample)
    setkey(rawDT, mz, sample)
    for (i in seq_along(irmzr)) {
        mzseq <- seq(start(irmzr[i]), end(irmzr[i]))
        for (s in .sampleNumber(object)) {
            rawDT[CJ(mzseq,s), shift := shiftsList[[i]][s], nomatch = 0]
            bgcorrDT[CJ(mzseq,s), shift := shiftsList[[i]][s], nomatch = 0]
        }
    }
    rawDT[, shift := ifelse(is.na(shift), 0, shift)]
    bgcorrDT[, shift := ifelse(is.na(shift), 0, shift)]
    rawDT[, scan := scan + shift]
    bgcorrDT[, scan := scan + shift]
    rawDT[, shift := NULL]
    bgcorrDT[, shift := NULL]
    rawDT <- rawDT[scan >= 1 & scan <= .maxScan(object)]
    bgcorrDT <- bgcorrDT[scan >= 1 & scan <= .maxScan(object)]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[rtAlignment] Remap scans .. done in %.1f secs.", stime))
    }
    .rawDT(object) <- rawDT
    .bgcorrDT(object) <- bgcorrDT
    .rtAlign(object) <- TRUE
    return(object)
}

densityEstimation <- function(object, dgridstep = dgridstep, dbandwidth = dbandwidth, 
                              outfileDens, verbose = FALSE) {
    stopifnot(is(object, "CMSproc"))
    bgcorrDT <- object@bgcorrDT
    getDensityEstimateApprox <- function(bgcorrDT, bw = dbandwidth,
                                         gridstep = dgridstep, maxbws = 4, mzParams) {
        gridseqMz <- seq(mzParams$minMZ, mzParams$maxMZ, gridstep[1])
        gridseqScan <- seq(1, mzParams$maxScan, gridstep[2])
        ## Assign each point to a grid location
        bgcorrDT[, gmz := cut(mz/1e5, breaks = gridseqMz, labels = FALSE)]
        bgcorrDT[, gscan := cut(scan, breaks = gridseqScan, labels = FALSE, include.lowest = TRUE)]
        ## Number of grid steps to meet maxbws (density beyond is approx 0)
        ng <- maxbws*bw/gridstep
        ## Sort by M/Z grid location than scan grid location
        setkey(bgcorrDT, gmz, gscan)
        tabgmz <- table(factor(bgcorrDT[,gmz], levels = seq_along(gridseqMz)))
        if(verbose) {
            message("[getDensityEstimateApprox] Getting sparse matrix entries (M/Z)")
        }
        ptime1 <- proc.time()
        spdensmz <- lapply(seq_along(gridseqMz), function(i) {
            whg <- max(i - ng[1], 1):min(i + ng[1], length(gridseqMz))
            if (i==1) {
                start <- 1
            } else {
                start <- sum(tabgmz[seq_len(whg[1]-1)])+1
            }
            end <- (start+sum(tabgmz[whg])-1)
            if (end < start)
                return(matrix(1,0,2))
            at <- start:end
            d <- rep(dnorm((gridseqMz[i]-gridseqMz[whg])/bw[1]), times = tabgmz[whg])
            return(cbind(at,d))
        })
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[getDensityEstimateApprox] Getting sparse matrix entries (M/Z) .. done in %.1f secs.", stime))
            message("[getDensityEstimateApprox] Constructing sparse matrix (M/Z)")
        }
        ptime1 <- proc.time()
        spmatmz <- sparseMatrix(
            i = do.call(c, lapply(seq_along(spdensmz), function(i) {
                               spdensmz[[i]][,1]
                           })),
            j = rep(seq_along(spdensmz), times = sapply(spdensmz, nrow)),
            x = do.call(c, lapply(seq_along(spdensmz), function(i) {
                               spdensmz[[i]][,2]
                           })),
            dims = c(nrow(bgcorrDT), length(spdensmz))
        )
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[getDensityEstimateApprox] Constructing sparse matrix (M/Z) .. done in %.1f secs.", stime))
        }
        rm(spdensmz)
        if(verbose) {
            message("[getDensityEstimateApprox] Getting sparse matrix entries (scan) + computing density")
        }
        gscan <- bgcorrDT[,gscan]
        intens <- bgcorrDT[,intensity]
        intens <- intens - min(intens)
        dtgscan <- data.table(gscan = gscan, row = seq_along(gscan))
        setkey(dtgscan, gscan)
        denom <- nrow(bgcorrDT)*prod(bw)*sum(bgcorrDT[,intensity])
        ptime1 <- proc.time()
        dens <- lapply(seq_along(gridseqScan), function(i) {
            whg <- max(i - ng[2], 1):min(i + ng[2], length(gridseqScan))
            at <- dtgscan[.(whg), row, nomatch = 0]
            ## Do intensity weighting here
            d <- dnorm((gridseqScan[i]-gscan[at])/bw[2])*bgcorrDT[at,intensity]
            ## rep(1) instead of rep(i) because we're forming one column matrices
            ## otherwise dims won't be correct
            spmatscan <- sparseMatrix(i = at, j = rep(1, length(at)), x = d, dims = c(nrow(spmatmz), 1))
            dens <- crossprod(spmatmz, spmatscan)/denom
            return(dens)
        })
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[getDensityEstimateApprox] Getting sparse matrix entries (scan) + computing density .. done in %.1f secs.", stime))
        }
        return(dens)
    }

    if(verbose) {
        message("[densityEstimation] Get density estimate")
    }
    ## First check to see if the density has been pre-computed
    if (!is.null(outfileDens)) {
        if (file.exists(outfileDens)) {
            if(verbose) {
                message("[getDensityEstimateApprox] Get density estimate .. loading")
            }
            load(outfileDens)
        } else {
            ptime1 <- proc.time()
            densList <- getDensityEstimateApprox(bgcorrDT = bgcorrDT, bw = dbandwidth,
                                                 gridstep = dgridstep, maxbws = 4, mzParams = object@mzParams)
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if(verbose) {
                message(sprintf("[densityEstimation] Get density estimate .. done in %.1f secs.", stime))
                message("[densityEstimation] Saving density estimate")
            }
            save(densList, file = outfileDens)
        }
    } else {
        ptime1 <- proc.time()
        densList <- getDensityEstimateApprox(bgcorrDT = bgcorrDT, bw = dbandwidth,
                                             gridstep = dgridstep, maxbws = 4, mzParams = object@mzParams)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[densityEstimation] Get density estimate  .. done in %.1f secs.", stime))
        }
    }
    
    if(verbose) {
        message("[densityEstimation] Make density matrix")
    }
    ptime1 <- proc.time()
    dmat <- do.call(cbind, lapply(densList, as.matrix))
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[densityEstimation] Make density matrix .. done in %.1f secs.", stime))
    }
    rm(densList)
    rownames(dmat) <- seq(.minMZ(object), .maxMZ(object), length.out = nrow(dmat))
    colnames(dmat) <- seq(1, .maxScan(object), dgridstep[2])
    return(list(dmat = dmat))
}

bakedpi <- function(cmsRaw, dbandwidth = c(0.005, 10),
                    dgridstep = c(0.005, 1), outfileDens = NULL,
                    dortalign = FALSE, mzsubset = NULL, verbose = TRUE) {
    ## Check that bandwidth is >= gridstep for binning purposes
    stopifnot(sum(dbandwidth >= dgridstep)==2)
    .isArgumentTwoVector(dbandwidth)
    .isArgumentTwoVector(dgridstep)

    ## Set verbosity options
    subverbose <- max(as.integer(verbose) - 1L, 0)

    ## Subset if specified
    cmsRaw <- .subsetByMZ(cmsRaw, mzsubset)
    
    if(verbose) {
        message("[bakedpi] Background correction")
    }
    object <- backgroundCorrection(object = cmsRaw, verbose = subverbose)

    if (dortalign) {
        object <- rtAlignment(object = object, verbose = subverbose)
    }

    if(verbose) {
        message("[bakedpi] Density estimation")
    }
    dmat <- densityEstimation(object = object, dbandwidth = dbandwidth,
                              dgridstep = dgridstep, outfileDens = outfileDens,
                              verbose = subverbose)$dmat
    .densityEstimate(object) <- dmat
    ## Compute and store density quantiles
    qs <- seq(0,1,0.001)
    .densityQuantiles(object) <- quantile(dmat[dmat!=0], qs)

    return(object)
}
