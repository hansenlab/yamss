
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

readRawDataAsDataTable <- function(obj, mzsubset = NULL, verbose = FALSE) {
    if(verbose) {
        message(sprintf("[readRawDataAsDataTable]: Reading %i files", length(obj@fileNames)))
    }
    obj@rawPeakInfo <- getPeakInfo(obj@fileNames)
    ## Make raw data matrix and data.table
    rawdatamat <- do.call(rbind, lapply(seq_along(obj@rawPeakInfo), function(s) {
                                     cbind(do.call(rbind, lapply(seq_along(obj@rawPeakInfo[[s]]), function(scan) {
                                                              cbind(obj@rawPeakInfo[[s]][[scan]], scan)
                                                          })), s)
                                 }))
    colnames(rawdatamat) <- c("mz", "intensity", "scan", "sample")
    DT <- data.table(mz = as.integer(rawdatamat[,"mz"]*1e5),
                     intensity = rawdatamat[,"intensity"],
                     scan = rawdatamat[,"scan"],
                     sample = rawdatamat[,"sample"])
    ## Subset if specified
    setkey(DT, mz, scan)
    if (!is.null(mzsubset)) {
        mzseq <- seq(as.integer(mzsubset[1]*1e5), as.integer(mzsubset[2]*1e5))
        DT <- DT[.(mzseq), nomatch = 0]
    }

    ## Get minimum and maximum M/Zs and scan numbers
    obj@mzParams <- list(maxScan = max(DT[,scan]),
                         maxMZ = 10*ceiling(max(DT[,mz])/1e5/10),
                         minMZ = 10*floor(min(DT[,mz])/1e5/10))
    obj@rawPeakDT <- DT
    return(obj)
}

backgroundCorrection <- function(obj, verbose = FALSE) {
    mzParams <- obj@mzParams
    DT <- obj@rawPeakDT
    setkey(DT, mz, scan, sample)
    mzbreaks <- c(seq(mzParams$minMZ, mzParams$maxMZ, 10), mzParams$maxMZ)
    scanbreaks <- seq(1, mzParams$maxScan, 40)
    scanbreaks[length(scanbreaks)] <- mzParams$maxScan
    if(verbose) {
        message("[backgroundCorrection] get marginal intensities")
    }
    ptime1 <- proc.time()
    densGrid <- lapply(1:(length(mzbreaks)-1), function(m) {
        mzseq <- seq(as.integer(mzbreaks[m]*1e5), as.integer(mzbreaks[m+1]*1e5))
        DTmz <- DT[.(mzseq), nomatch = 0]
        lapply(1:(length(scanbreaks)-1), function(rt) {
            DTmzscan <- DTmz[scan %in% scanbreaks[rt]:scanbreaks[rt+1]]
            lapply(seq_along(obj@fileNames), function(s) {
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
        message(sprintf("[backgroundCorrection]   .. done in %.1f secs.", stime))
    }
    ## Estimate retention time window-specific background levels
    bgsd <- 0 # SD of normal distribution characterizing noise intensities
    r <- dnorm(bgsd)/dnorm(0)
    bgmeans <- lapply(1:(length(mzbreaks)-1), function(m) {
        lapply(1:(length(scanbreaks)-1), function(rt) {
            lapply(seq_along(obj@fileNames), function(s) {
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
        message("[backgroundCorrection] get region-specific background trends")
    }
    ptime1 <- proc.time()
    smooths <- lapply(seq_along(obj@fileNames), function(s) {
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
        message(sprintf("[backgroundCorrection]  .. done in %.1f secs.", stime))
    }
    obj@bgSmooths <- smooths

    ## Perform background correction
    setkey(DT, mz, sample)
    if(verbose) {
        message("[backgroundCorrection] correct intensities")
    }
    ptime1 <- proc.time()
    DTbgcorr <- rbindlist(lapply(seq_along(obj@fileNames), function(s) {
        bgmeanmatSmoothed <- smooths[[s]]
        mzbounds <- attr(bgmeanmatSmoothed, "mzbounds")
        rbindlist(lapply(1:nrow(mzbounds), function(i) {
            mzseq <- seq(as.integer(mzbounds[i,1]*1e5), as.integer(mzbounds[i,2]*1e5))
            bgtrend <- bgmeanmatSmoothed[,i]
            bgtrend[is.na(bgtrend)] <- 0
            DT[.(mzseq), nomatch = 0][sample==s][,bg := (2^bgtrend[scan])-1]
        }))
    }))
    DTbgcorr[, intensity := intensity - bg]
    DTbgcorr[, bg := NULL]
    DTbgcorr <- DTbgcorr[intensity > 0]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[backgroundCorrection]  .. done in %.1f secs.", stime))
    }
    obj@bgcorrDT <- DTbgcorr
    return(obj)
}

rtAlignment <- function(obj, verbose = FALSE) {
    mzParams <- obj@mzParams
    DT <- obj@rawPeakDT
    DTbgcorr <- obj@bgcorrDT
    if(verbose) {
        message("[rtAlignment] Get rough M/Z regions to align")
    }
    ## Get density of M/Z values
    ## This density will be thresholded to yield M/Z regions
    mzdens <- density(DTbgcorr[,mz]/1e5, weights = DTbgcorr[,intensity]/sum(DTbgcorr[,intensity]), n = 2^ceiling(log2(nrow(DTbgcorr))), bw = 0.005)
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
        message("[rtAlignment] Get XICs for these regions")
    }
    ptime1 <- proc.time()
    xics <- getAllXics(mzranges = irmzr, DT = DT)
    scans <- 1:obj@mzParams$maxScan
    xicsImputed <- lapply(xics, function(x) {
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
        message(sprintf("[rtAlignment]  .. done in %.1f secs.", stime))
        message("[rtAlignment] Find best shifts")
    }
    ptime1 <- proc.time()
    shifts <- -20:20
    DTbgcorr[, scanorig := scan]
    DT[, scanorig := scan]
    shiftsList <- lapply(seq_along(irmzr), function(i) {
        xicmat <- xics[[i]]
        xicimpmat <- xicsImputed[[i]]
        refsamp <- which.max(colSums(xicmat))
        bestShiftBySample <- sapply(seq_along(obj@fileNames), function(s) {
            if (s==refsamp) {
                return(0)
            }
            corrShifts <- sapply(shifts, function(shift) {
                if (shift < 0) {
                    x <- tail(xicimpmat[,s], shift)
                    ref <- head(xicimpmat[,refsamp], shift)
                } else if (shift==0) {
                    x <- xicimpmat[,s]
                    ref <- xicimpmat[,refsamp]
                } else {
                    x <- head(xicimpmat[,s], -shift)
                    ref <- tail(xicimpmat[,refsamp], -shift)
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
        message(sprintf("[rtAlignment]  .. done in %.1f secs.", stime))
        message("[rtAlignment] Remap scans")
    }
    obj@alignments <- shiftsList
    ptime1 <- proc.time()
    setkey(DTbgcorr, mz, sample)
    setkey(DT, mz, sample)
    for (i in seq_along(irmzr)) {
        mzseq <- seq(start(irmzr[i]), end(irmzr[i]))
        for (s in seq_along(obj@fileNames)) {
            DT[CJ(mzseq,s), shift := shiftsList[[i]][s], nomatch = 0]
            DTbgcorr[CJ(mzseq,s), shift := shiftsList[[i]][s], nomatch = 0]
        }
    }
    DT[, shift := ifelse(is.na(shift), 0, shift)]
    DTbgcorr[, shift := ifelse(is.na(shift), 0, shift)]
    DT[, scan := scan + shift]
    DTbgcorr[, scan := scan + shift]
    DT[, shift := NULL]
    DTbgcorr[, shift := NULL]
    DT <- DT[scan >= 1 & scan <= obj@mzParams$maxScan]
    DTbgcorr <- DTbgcorr[scan >= 1 & scan <= obj@mzParams$maxScan]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf(".. done in %.1f secs.", stime))
    }
    obj@rawPeakDT <- DT
    obj@bgcorrDT <- DTbgcorr
    return(obj)
}

densityEstimation <- function(obj, dgridstep = dgridstep, dbandwidth = dbandwidth, 
                              outfileDens, verbose = FALSE) {
    DTbgcorr <- obj@bgcorrDT
    getDensityEstimateApprox <- function(DTbgcorr, bw = dbandwidth,
                                         gridstep = dgridstep, maxbws = 4, mzParams) {
        gridseqMz <- seq(mzParams$minMZ, mzParams$maxMZ, gridstep[1])
        gridseqScan <- seq(1, mzParams$maxScan, gridstep[2])
        ## Assign each point to a grid location
        DTbgcorr[, gmz := cut(mz/1e5, breaks = gridseqMz, labels = FALSE)]
        DTbgcorr[, gscan := cut(scan, breaks = gridseqScan, labels = FALSE, include.lowest = TRUE)]
        ## Number of grid steps to meet maxbws (density beyond is approx 0)
        ng <- maxbws*bw/gridstep
        ## Sort by M/Z grid location than scan grid location
        setkey(DTbgcorr, gmz, gscan)
        tabgmz <- table(factor(DTbgcorr[,gmz], levels = 1:length(gridseqMz)))
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
            message(sprintf("[getDensityEstimateApprox]  .. done in %.1f secs.", stime))
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
            dims = c(nrow(DTbgcorr), length(spdensmz))
        )
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[getDensityEstimateApprox]  .. done in %.1f secs.", stime))
        }
        rm(spdensmz)
        if(verbose) {
            message("[getDensityEstimateApprox] Getting sparse matrix entries (scan) + computing density")
        }
        gscan <- DTbgcorr[,gscan]
        intens <- DTbgcorr[,intensity]
        intens <- intens - min(intens)
        dtgscan <- data.table(gscan = gscan, row = seq_along(gscan))
        setkey(dtgscan, gscan)
        denom <- nrow(DTbgcorr)*prod(bw)*sum(DTbgcorr[,intensity])
        ptime1 <- proc.time()
        dens <- lapply(seq_along(gridseqScan), function(i) {
            whg <- max(i - ng[2], 1):min(i + ng[2], length(gridseqScan))
            at <- dtgscan[.(whg), row, nomatch = 0]
            ## Do intensity weighting here
            d <- dnorm((gridseqScan[i]-gscan[at])/bw[2])*DTbgcorr[at,intensity]
            ## rep(1) instead of rep(i) because we're forming one column matrices
            ## otherwise dims won't be correct
            spmatscan <- sparseMatrix(i = at, j = rep(1, length(at)), x = d, dims = c(nrow(spmatmz), 1))
            dens <- crossprod(spmatmz, spmatscan)/denom
            return(dens)
        })
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[getDensityEstimateApprox]  .. done in %.1f secs.", stime))
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
                message("[getDensityEstimateApprox]  .. loading")
            }
            load(outfileDens)
        } else {
            ptime1 <- proc.time()
            densList <- getDensityEstimateApprox(DTbgcorr = DTbgcorr, bw = dbandwidth,
                                                 gridstep = dgridstep, maxbws = 4, mzParams = obj@mzParams)
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if(verbose) {
                message(sprintf("[densityEstimation] estimation .. done in %.1f secs.", stime))
                message("[densityEstimation] Saving density estimate")
            }
            save(densList, file = outfileDens)
        }
    } else {
        ptime1 <- proc.time()
        densList <- getDensityEstimateApprox(DTbgcorr = DTbgcorr, bw = dbandwidth,
                                             gridstep = dgridstep, maxbws = 4, mzParams = obj@mzParams)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) {
            message(sprintf("[densityEstimation] estimation .. done in %.1f secs.", stime))
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
        message(sprintf("[densityEstimation] .. done in %.1f secs.", stime))
    }
    rm(densList)
    rownames(dmat) <- seq(obj@mzParams$minMZ, obj@mzParams$maxMZ, length.out = nrow(dmat))
    colnames(dmat) <- seq(1, obj@mzParams$maxScan, dgridstep[2])
    return(list(dmat = dmat))
}

getCutoff <- function(obj, densmat, by = 2, verbose = FALSE) {
    if(verbose) {
        message("[getDensityCutoff] Get density cutoff")
    }
    ptime1 <- proc.time()
    mzParams <- obj@mzParams
    DTbgcorr <- obj@bgcorrDT
    mzregions <- seq(mzParams$minMZ, mzParams$maxMZ, by)
    setkey(DTbgcorr, mz)
    densmatMzs <- as.numeric(rownames(densmat))
    densmatScans <- as.numeric(colnames(densmat))
    ## For each M/Z region, find the data point with the highest intensity
    ## Find the peak region corrresponding to this data point
    ## Get the density estimate values in this region
    dlist <- lapply(1:(length(mzregions)-1), function(i) {
        ## Get the point with the highest intensity in this region
        mzseq <- seq(mzregions[i]*1e5, mzregions[i+1]*1e5-1)
        subdt <- DTbgcorr[.(mzseq), nomatch = 0]
        if (nrow(subdt)==0)
            return(c())
        whm <- which.max(subdt[,intensity])
        ## Define a small region about this point in which to roughly define a peak area
        scanwindow <- c(floor(subdt[whm,scan]-5), ceiling(subdt[whm,scan]+5))
        scan <- subdt[whm,scan]
        mz <- subdt[whm,mz]/1e5
        mzwindow <- c(mz-(mz*30/1e6), mz+(mz*30/1e6))
        mzseq <- seq(as.integer(mzwindow[1]*1e5), as.integer(mzwindow[2]*1e5))
        subdt <- DTbgcorr[.(mzseq), nomatch = 0][scan >= scanwindow[1] & scan <= scanwindow[2]]
        ## Skip if number of points in region isn't >= number of samples
        if (nrow(subdt) < length(obj@fileNames))
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
        subdt <- DTbgcorr[.(mzseq), nomatch = 0]
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
        message(sprintf("[getDensityCutoff]  .. done in %.1f secs.", stime))
    }
    return(cutoff)
}

getBlobs <- function(densmat, dcutoff, verbose = FALSE) {
    if(verbose) {
        message("[getBlobs] Getting blobs")
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
    dtblobs <- data.table(row = wh[,1], col = scans[wh[,2]], blobnum = blobs[wh])
    setkey(dtblobs, blobnum)
    blobsDT <- dtblobs[, .(mzmin = mzs[min(row)], mzmax = mzs[max(row)+1], scanmin = min(col), scanmax = max(col)), by = blobnum]
    ## Make matrix of blob information
    blobs <- cbind(blobsDT[,mzmin], blobsDT[,mzmax], blobsDT[,scanmin], blobsDT[,scanmax], blobsDT[,blobnum])
    colnames(blobs) <- c("mzmin", "mzmax", "scan.min", "scan.max", "blobnum")
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[getBlobs]  .. done in %.1f secs.", stime))
    }
    return(blobs)
}

getXICsAndQuantifyWithRetentionTime <- function(obj, verbose = FALSE) {
    if(verbose) {
        message("[getXICsAndQuantifyWithRetentionTime] compute XICs")
    }
    mzParams <- obj@mzParams
    DT <- obj@rawPeakDT
    setkey(DT, mz, scan)
    ptime1 <- proc.time()
    obj@xicsRaw <- lapply(1:nrow(obj@peakBounds), function(i) {
        mzseq <- seq(as.integer(obj@peakBounds[i,"mzmin"]*1e5), as.integer(obj@peakBounds[i,"mzmax"]*1e5))
        dt <- DT[.(mzseq), nomatch = 0]
        dt <- dt[, xic := log2(max(intensity)+1), by = .(scan, sample)]
        dup <- duplicated(dt[,.(scan,sample)])
        dt <- dt[!dup]
        xics <- lapply(seq_along(obj@fileNames), function(s) {
            subdt <- dt[sample==s]
            if (nrow(subdt) < 2) {
                return(approxfun(x = 1:2, y = rep(0,2), rule = 2))
            } else {
                return(approxfun(x = subdt[,scan], y = subdt[,xic], rule = 2))
            }
        })
        return(xics)
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
    quantmat <- do.call(rbind, lapply(seq_along(obj@xicsRaw), function(i) {
                                   wh <- which.min(abs(scanseq-obj@peakBounds[i,"scan.min"])):which.min(abs(scanseq-obj@peakBounds[i,"scan.max"]))
                                   sapply(obj@xicsRaw[[i]], function(xic) {
                                       f <- (2^xic(scanseq))-1
                                       return(sum(f[wh])*scanstep)
                                   })
                               }))
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[getXICsAndQuantifyWithRetentionTime]  .. done in %.1f secs.", stime))
    }
    obj@peakQuants <- quantmat
    return(obj)
}

getXICsAndQuantifyWithoutRetentionTime <- function(obj, verbose = FALSE) {
    if(verbose) {
        message("[getXICsAndQuantifyWithoutRetentionTime] compute XICs")
    }
    mzParams <- obj@mzParams
    DT <- obj@rawPeakDT
    ptime1 <- proc.time()
    obj@xicsRaw <- getAllXics(mzranges = IRanges(start = as.integer(obj@peakBounds[,"mzmin"]*1e5),
                                                 end = as.integer(obj@peakBounds[,"mzmax"]*1e5)), DT = DT)
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
    obj@xicsImputed <- lapply(obj@xicsRaw, function(x) {
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
    quantmat <- do.call(rbind, lapply(seq_along(obj@xicsImputed), function(i) {
                                   mat <- obj@xicsImputed[[i]][obj@peakBounds[i, "scan.min"]:min(c(obj@peakBounds[i, "scan.max"], mzParams$maxScan)),,drop = FALSE]
                                   mat <- (2^mat)-1
                                   colSums(mat)
                               }))
    obj@peakQuants <- quantmat
    return(obj)
}

bakedpi <- function(files, dbandwidth = c(0.005, 10),
                    dgridstep = c(0.005, 1), outfileDens = NULL,
                    dortalign = FALSE, mzsubset = NULL, verbose = TRUE) {
    ## Check that bandwidth is >= gridstep for binning purposes
    stopifnot(sum(dbandwidth >= dgridstep)==2)
    .isArgumentTwoVector(dbandwidth)
    .isArgumentTwoVector(dgridstep)
    subverbose <- max(as.integer(verbose) - 1L, 0)
    obj <- new("CMS", fileNames = files, rtalign = dortalign)
    ## Parse raw data
    if(verbose) {
        message("[bakedpi] Reading data")
    }
    obj <- readRawDataAsDataTable(obj = obj, mzsubset = mzsubset, verbose = subverbose)

    if(verbose) {
        message("[bakedpi] Background correction")
    }
    obj <- backgroundCorrection(obj = obj, verbose = subverbose)

    if (dortalign) {
        obj <- rtAlignment(obj = obj, verbose = subverbose)
    }

    if(verbose) {
        message("[bakedpi] Density estimation")
    }
    dmat <- densityEstimation(obj = obj, dbandwidth = dbandwidth,
                              dgridstep = dgridstep, outfileDens = outfileDens,
                              verbose = subverbose)$dmat
    obj@density <- dmat

    if(verbose) {
        message("[bakedpi] Computing cutoff")
    }
    cutoff <- getCutoff(obj = obj, densmat = dmat, by = 2, verbose = subverbose)
    qs <- seq(0,1,0.001)
    obj@densityQuants <- quantile(dmat[dmat!=0], qs)
    qcutoff <- which.min(abs(cutoff-obj@densityQuants))
    qref <- which.min(abs(qs-0.99))
    obj@densityCutoff <- max(cutoff, obj@densityQuants[qref])
    if(verbose) {
        message("[bakedpi] Getting blobs")
    }
    obj@peakBounds <- getBlobs(densmat = dmat, dcutoff = obj@densityCutoff, verbose = subverbose)

    ## Get XICs and quantifications - methods differ if RT alignment was performed
    if(verbose) {
        message("[bakedpi] Quantifying")
    }
    if (dortalign) {
        obj <- getXICsAndQuantifyWithRetentionTime(obj = obj, verbose = subverbose)
    } else {
        obj <- getXICsAndQuantifyWithoutRetentionTime(obj = obj, verbose = subverbose)
    }
    
    ## Differential analysis
    ## message("[bakedpi] Differential analysis")
    ## diffrep <- getDiffTable(obj@peakQuants, classes)
    ## obj@diffrep <- diffrep
    return(obj)
}
