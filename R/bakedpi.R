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


.setMZParams <- function(rawDT) {
    mzParams <- list(
        mminScan = min(rawDT[,scan]),
        maxScan = max(rawDT[,scan]),
        minMZraw = min(rawDT[,mz])/1e5,
        maxMZraw = max(rawDT[,mz])/1e5,
        minMZ = 10*floor(min(rawDT[,mz])/1e6),
        maxMZ = 10*ceiling(max(rawDT[,mz])/1e6))
    mzParams
}
    
.subsetByMZ <- function(object, mzsubset = NULL) {
    if(is.null(mzsubset))
        return(object)
    rawDT <- object@rawDT
    setkey(rawDT, mz, scan)
    mzseq <- seq(as.integer(mzsubset[1]*1e5), as.integer(mzsubset[2]*1e5))
    object@rawDT <- rawDT[.(mzseq), nomatch = 0]
    object@mzParams <- .setMZParams(object@rawDT)
    object
}


readMSdata <- function(files, colData = NULL, mzsubset = NULL, verbose = FALSE) {
    if(verbose) {
        message(sprintf("[readRaw]: Reading %i files", length(files)))
    }
    ## Check that file exists
    stopifnot(all(file.exists(files)))
    ## If colData is supplied, check nrow(colData)==length(files)
    if (!is.null(colData)) {
        stopifnot(nrow(colData)==length(files))
    }
    cmsRaw <- new("CMSraw")
    rawPeakInfo <- getPeakInfo(files)
    ## Make raw data matrix and data.table
    rawdatamat <- do.call(rbind, lapply(seq_along(rawPeakInfo), function(s) {
                                     cbind(do.call(rbind, lapply(seq_along(rawPeakInfo[[s]]), function(scan) {
                                                              cbind(rawPeakInfo[[s]][[scan]], scan)
                                                          })), s)
                                 }))
    colnames(rawdatamat) <- c("mz", "intensity", "scan", "sample")
    rawDT <- data.table(mz = as.integer(rawdatamat[,"mz"]*1e5),
                     intensity = rawdatamat[,"intensity"],
                     scan = rawdatamat[,"scan"],
                     sample = rawdatamat[,"sample"])
    ## Get minimum and maximum M/Zs and scan numbers
    cmsRaw@mzParams <- .setMZParams(rawDT)
    cmsRaw@rawDT <- rawDT
    ## Make coltype data frame
    ## FIXME: check that colData, if not null, does not have columns names "files" and "sample"
    fileData <- DataFrame(sample = seq_along(files), files = files)
    if(is.null(colData)) {
        cmsRaw@colData <- fileData
    } else {
        cmsRaw@colData <- cbind(fileData, colData)
    }
    cmsRaw <- .subsetByMZ(cmsRaw, mzsubset = mzsubset)
    return(cmsRaw)
}

backgroundCorrection <- function(object, verbose = FALSE) {
    stopifnot(is(object, "CMSraw"))
    mzParams <- object@mzParams
    rawDT <- object@rawDT
    setkey(rawDT, mz, scan, sample)
    mzbreaks <- c(seq(mzParams$minMZ, mzParams$maxMZ, 10), mzParams$maxMZ)
    scanbreaks <- seq(1, mzParams$maxScan, 40)
    scanbreaks[length(scanbreaks)] <- mzParams$maxScan
    if(verbose) {
        message("[backgroundCorrection] Get marginal intensities")
    }
    ptime1 <- proc.time()
    densGrid <- lapply(1:(length(mzbreaks)-1), function(m) {
        mzseq <- seq(as.integer(mzbreaks[m]*1e5), as.integer(mzbreaks[m+1]*1e5))
        DTmz <- rawDT[.(mzseq), nomatch = 0]
        lapply(1:(length(scanbreaks)-1), function(rt) {
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
    bgmeans <- lapply(1:(length(mzbreaks)-1), function(m) {
        lapply(1:(length(scanbreaks)-1), function(rt) {
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
        bgmeanmatThisSample <- bgmeanmatThisSample[,keepcols]
        ## Each col is a smoothed background trend across scans for a particular M/Z region
        bgmeanmatSmoothed <- do.call(cbind, lapply(1:ncol(bgmeanmatThisSample), function(i) {
                                                dists <- (1:ncol(bgmeanmatThisSample))-i
                                                wts <- dnorm(dists/4)
                                                weightmat <- matrix(wts, nrow = nrow(bgmeanmatThisSample), ncol = ncol(bgmeanmatThisSample), byrow = TRUE)
                                                df <- data.frame(intens = as.numeric(bgmeanmatThisSample), scan = rep(head(scanbreaks, -1), times = ncol(bgmeanmatThisSample)), weight = as.numeric(weightmat))
                                                df <- df[complete.cases(df),]
                                                lofit <- loess(intens ~ scan, data = df, weights = weight, span = 0.1)
                                                predict(lofit, 1:mzParams$maxScan)
                                            }))
        mzbounds <- cbind(head(mzbreaks, -1), tail(mzbreaks, -1))
        mzbounds <- mzbounds[keepcols,]
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
        rbindlist(lapply(1:nrow(mzbounds), function(i) {
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
    out@bgcorrDT <- bgcorrDT
    return(out)
}

rtAlignment <- function(object, verbose = FALSE) {
    stopifnot(is(object, "CMSproc"))
    mzParams <- object@mzParams
    rawDT <- object@rawDT
    bgcorrDT <- object@bgcorrDT
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
    scans <- 1:object@mzParams$maxScan
    eicsImputed <- lapply(eics, function(x) {
        do.call(cbind, lapply(1:ncol(x), function(col) {
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
    rawDT <- rawDT[scan >= 1 & scan <= object@mzParams$maxScan]
    bgcorrDT <- bgcorrDT[scan >= 1 & scan <= object@mzParams$maxScan]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[rtAlignment] Remap scans .. done in %.1f secs.", stime))
    }
    object@rawDT <- rawDT
    object@bgcorrDT <- bgcorrDT
    object@rtAlign <- TRUE
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
        tabgmz <- table(factor(bgcorrDT[,gmz], levels = 1:length(gridseqMz)))
        if(verbose) {
            message("[getDensityEstimateApprox] Getting sparse matrix entries (M/Z)")
        }
        ptime1 <- proc.time()
        spdensmz <- lapply(seq_along(gridseqMz), function(i) {
            whg <- max(i - ng[1], 1):min(i + ng[1], length(gridseqMz))
            if (i==1) {
                start <- 1
            } else {
                start <- sum(tabgmz[1:(whg[1]-1)])+1
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
    rownames(dmat) <- seq(object@mzParams$minMZ, object@mzParams$maxMZ, length.out = nrow(dmat))
    colnames(dmat) <- seq(1, object@mzParams$maxScan, dgridstep[2])
    return(list(dmat = dmat))
}

getCutoff <- function(object, by = 2, verbose = FALSE) {
    if(verbose) {
        message("[getDensityCutoff] Get density cutoff")
    }
    stopifnot(is(object, "CMSproc"))
    ptime1 <- proc.time()
    mzParams <- object@mzParams
    bgcorrDT <- object@bgcorrDT
    densmat <- object@density
    mzregions <- seq(mzParams$minMZ, mzParams$maxMZ, by)
    setkey(bgcorrDT, mz)
    densmatMzs <- as.numeric(rownames(densmat))
    densmatScans <- as.numeric(colnames(densmat))
    ## For each M/Z region, find the data point with the highest intensity
    ## Find the peak region corrresponding to this data point
    ## Get the density estimate values in this region
    dlist <- lapply(1:(length(mzregions)-1), function(i) {
        ## Get the point with the highest intensity in this region
        mzseq <- seq(mzregions[i]*1e5, mzregions[i+1]*1e5-1)
        subdt <- bgcorrDT[.(mzseq), nomatch = 0]
        if (nrow(subdt)==0)
            return(c())
        whm <- which.max(subdt[,intensity])
        ## Define a small region about this point in which to roughly define a peak area
        scanwindow <- c(floor(subdt[whm,scan]-5), ceiling(subdt[whm,scan]+5))
        scan <- subdt[whm,scan]
        mz <- subdt[whm,mz]/1e5
        mzwindow <- c(mz-(mz*30/1e6), mz+(mz*30/1e6))
        mzseq <- seq(as.integer(mzwindow[1]*1e5), as.integer(mzwindow[2]*1e5))
        subdt <- bgcorrDT[.(mzseq), nomatch = 0][scan >= scanwindow[1] & scan <= scanwindow[2]]
        ## Skip if number of points in region isn't >= number of samples
        if (nrow(subdt) < length(.sampleNumber(object)))
            return(c())
        ## Get the M/Z range of this rough "peak"
        dens <- density(subdt[,mz]/1e5)
        maxIndex <- which.max(dens$y)
        left <- maxIndex - which.max(dens$y[maxIndex:1] < max(dens$y)*0.25) + 1
        right <- maxIndex + which.max(dens$y[maxIndex:length(dens$y)] < max(dens$y)*0.25) - 1
        mzwindow <- dens$x[c(left, right)]
        densmatRows <- which.min(abs(densmatMzs-mzwindow[1])):which.min(abs(densmatMzs-mzwindow[2]))
        ## Get the scan range of this rough "peak" (left, right)
        mzseq <- seq(as.integer(mzwindow[1]*1e5), as.integer(mzwindow[2]*1e5))
        subdt <- bgcorrDT[.(mzseq), nomatch = 0]
        subdt <- subdt[, .N, by = scan]
        numperscan <- rep(0, mzParams$maxScan)
        numperscan[subdt[,scan]] <- subdt[,N]
        names(numperscan) <- 1:mzParams$maxScan
        numperscan <- numperscan[as.character(densmatScans)]
        scanIndex <- which.min(abs(densmatScans-scan))
        left <- ifelse(sum(numperscan[scanIndex:1]==0)==0, 1, scanIndex-which.max(numperscan[scanIndex:1]==0)+1)
        right <- ifelse(sum(numperscan[scanIndex:length(densmatScans)]==0)==0,
                        length(densmatScans),
                        scanIndex+which.max(numperscan[scanIndex:length(densmatScans)]==0)-1)
        return(as.numeric(densmat[densmatRows,left:right]))
    })
    ## Use features of the distribution of these density values to select a cutoff
    ## Obtain quantiles of these density values in the different regions
    ## Look at the modes of these distributions - for most quantiles, the mode is near 0
    ## Choose the cutoff as the first mode that "jumps" away from zero
    qs <- seq(0.001,0.999,0.001)
    modes <- sapply(qs, function(qu) {
        quants <- sapply(dlist, quantile, qu)
        d <- density(quants, from = min(quants, na.rm = TRUE), na.rm = TRUE)
        d$x[which.max(d$y)]
    })
    sortedModes <- sort(modes)
    r <- tail(sortedModes, -1)/head(sortedModes, -1)
    cp <- which.max(r > 10)+1
    cutoff <- sortedModes[cp]

    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[getDensityCutoff] Get density cutoff .. done in %.1f secs.", stime))
    }
    return(cutoff)
}

computePeakBounds <- function(densmat, dcutoff, verbose = FALSE) {
    if(verbose) {
        message("[computePeakBounds] Computing peak bounds")
    }
    ptime1 <- proc.time()
    ## Get grid ticks
    mzs <- as.numeric(rownames(densmat))
    mzs <- c(mzs, tail(mzs, 1)+(mzs[2]-mzs[1]))
    scans <- as.numeric(colnames(densmat))
    ## Connected-components labeling
    bool <- densmat > dcutoff
    blobs <- bwlabel(bool)
    wh <- which(blobs!=0, arr.ind = TRUE)
    dtblobs <- data.table(row = wh[,1], col = scans[wh[,2]], peaknum = blobs[wh])
    setkey(dtblobs, peaknum)
    blobsDT <- dtblobs[, .(mzmin = mzs[min(row)], mzmax = mzs[max(row)+1], scanmin = min(col), scanmax = max(col)), by = peaknum]
    ## Make matrix of blob information
    blobs <- cbind(blobsDT[,mzmin], blobsDT[,mzmax], blobsDT[,scanmin], blobsDT[,scanmax], blobsDT[,peaknum])
    colnames(blobs) <- c("mzmin", "mzmax", "scan.min", "scan.max", "peaknum")
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[computePeakBounds] Computing peak bounds .. done in %.1f secs.", stime))
    }
    return(blobs)
}

getXICsAndQuantifyWithRetentionTime <- function(object, verbose = FALSE) {
    if(verbose) {
        message("[getXICsAndQuantifyWithRetentionTime] compute XICs")
    }
    stopifnot(is(object, "CMSproc"))
    mzParams <- object@mzParams
    rawDT <- object@rawDT
    setkey(rawDT, mz, scan)
    ptime1 <- proc.time()
    eicsRaw <- lapply(1:nrow(object@peakBounds), function(i) {
        mzseq <- seq(as.integer(object@peakBounds[i,"mzmin"]*1e5), as.integer(object@peakBounds[i,"mzmax"]*1e5))
        dt <- rawDT[.(mzseq), nomatch = 0]
        dt <- dt[, eic := log2(max(intensity)+1), by = .(scan, sample)]
        dup <- duplicated(dt[,.(scan,sample)])
        dt <- dt[!dup]
        eics <- lapply(.sampleNumber(object), function(s) {
            subdt <- dt[sample==s]
            if (nrow(subdt) < 2) {
                return(approxfun(x = 1:2, y = rep(0,2), rule = 2))
            } else {
                return(approxfun(x = subdt[,scan], y = subdt[,eic], rule = 2))
            }
        })
        return(eics)
    })
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf(".. done in %.1f secs.", stime))
        message("[getXICsAndQuantifyWithRetentionTime] quantify")
    }
    scanstep <- 0.01
    scanseq <- seq(0, mzParams$maxScan, scanstep)
    ptime1 <- proc.time()
    quantmat <- do.call(rbind, lapply(seq_along(eicsRaw), function(i) {
                                   wh <- which.min(abs(scanseq-object@peakBounds[i,"scan.min"])):which.min(abs(scanseq-object@peakBounds[i,"scan.max"]))
                                   sapply(eicsRaw[[i]], function(eic) {
                                       f <- (2^eic(scanseq))-1
                                       return(sum(f[wh])*scanstep)
                                   })
                               }))
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[getXICsAndQuantifyWithRetentionTime]  .. done in %.1f secs.", stime))
    }
    object@peakQuants <- quantmat
    return(object)
}

getXICsAndQuantifyWithoutRetentionTime <- function(object, verbose = FALSE) {
    if(verbose) {
        message("[getXICsAndQuantifyWithoutRetentionTime] compute XICs")
    }
    stopifnot(is(object, "CMSproc"))
    mzParams <- object@mzParams
    ptime1 <- proc.time()
    irmzr <- IRanges(start = as.integer(object@peakBounds[,"mzmin"]*1e5), 
                     end = as.integer(object@peakBounds[,"mzmax"]*1e5))
    eicsRaw <- getEICS(object, mzranges = irmzr)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf(".. done in %.1f secs.", stime))
    }
    scans <- 1:mzParams$maxScan

    if(verbose) {
        message("[getXICsAndQuantifyWithoutRetentionTime] Impute")
    }
    ptime1 <- proc.time()
    eicsImputed <- lapply(eicsRaw, function(x) {
        do.call(cbind, lapply(1:ncol(x), function(col) {
                           bool <- x[,col] > 1e-6
                           if (sum(bool) < 2)
                               return(x[,col])
                           approx(scans[bool], x[bool,col], xout = scans, rule = 2)$y
                       }))
    })
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[getXICsAndQuantifyWithoutRetentionTime]  .. done in %.1f secs.", stime))
    }
    
    if(verbose) {
        message("[getXICsAndQuantifyWithoutRetentionTime] quantify")
    }
    quantmat <- do.call(rbind, lapply(seq_along(eicsImputed), function(i) {
                                   mat <- eicsImputed[[i]][object@peakBounds[i, "scan.min"]:min(c(object@peakBounds[i, "scan.max"], mzParams$maxScan)),,drop = FALSE]
                                   mat <- (2^mat)-1
                                   colSums(mat)
                               }))
    object@peakQuants <- quantmat
    return(object)
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
    ## Fixme: let this be a subsetting function, which we may want to export
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
    object@density <- dmat

    if(verbose) {
        message("[bakedpi] Computing cutoff")
    }
    cutoff <- getCutoff(object = object, by = 2, verbose = subverbose)
    qs <- seq(0,1,0.001)
    object@densityQuantiles <- quantile(dmat[dmat!=0], qs)
    qcutoff <- which.min(abs(cutoff-object@densityQuantiles))
    qref <- which.min(abs(qs-0.99))
    object@densityCutoff <- max(cutoff, object@densityQuantiles[qref])
    if(verbose) {
        message("[bakedpi] Computing peak bounds")
    }
    object@peakBounds <- computePeakBounds(densmat = dmat, dcutoff = object@densityCutoff, verbose = subverbose)

    ## Get XICs and quantifications - methods differ if RT alignment was performed
    if(verbose) {
        message("[bakedpi] Quantifying peaks")
    }
    if (dortalign) {
        object <- getXICsAndQuantifyWithRetentionTime(object = object, verbose = subverbose)
    } else {
        object <- getXICsAndQuantifyWithoutRetentionTime(object = object, verbose = subverbose)
    }
    
    return(object)
}
