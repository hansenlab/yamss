getCutoff <- function(object, mzSpacing = 2, verbose = FALSE) {
    if(verbose) {
        message("[getDensityCutoff] Get density cutoff")
    }
    stopifnot(is(object, "CMSproc"))
    ptime1 <- proc.time()
    bgcorrDT <- .bgcorrDT(object)
    densmat <- densityEstimate(object)
    mzregions <- seq(.minMZ(object), .maxMZ(object), by = mzSpacing)
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
        numperscan <- rep(0, .maxScan(object))
        numperscan[subdt[,scan]] <- subdt[,N]
        names(numperscan) <- 1:.maxScan(object)
        numperscan <- numperscan[as.character(densmatScans)]
        scanIndex <- which.min(abs(densmatScans-scan))
        if (sum(numperscan[scanIndex:1]==0)==0) {
            left <- 1
        } else {
            left <- scanIndex-which.max(numperscan[scanIndex:1]==0)+1
        }
        if (sum(numperscan[scanIndex:length(densmatScans)]==0)==0) {
            right <- length(densmatScans)
        } else {
            right <- scanIndex+which.max(numperscan[scanIndex:length(densmatScans)]==0)-1
        }
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
    colnames(blobs) <- c("mzmin", "mzmax", "scanmin", "scanmax", "peaknum")
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[computePeakBounds] Computing peak bounds .. done in %.1f secs.", stime))
    }
    return(blobs)
}

getEICsAndQuantify <- function(object, peakBounds, verbose = FALSE) {
    if(verbose) {
        message("[getEICsAndQuantify] compute EICs")
    }
    stopifnot(is(object, "CMSproc"))
    rawDT <- .rawDT(object)
    setkey(rawDT, mz, scan)
    ptime1 <- proc.time()
    eicsRaw <- lapply(1:nrow(peakBounds), function(i) {
        mzseq <- seq(as.integer(peakBounds[i,"mzmin"]*1e5), as.integer(peakBounds[i,"mzmax"]*1e5))
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
        message("[getEICsAndQuantify] quantify")
    }
    scanstep <- 0.01
    scanseq <- seq(0, .maxScan(object), scanstep)
    ptime1 <- proc.time()
    quantmat <- do.call(rbind, lapply(seq_along(eicsRaw), function(i) {
                                   wh <- which.min(abs(scanseq-peakBounds[i,"scanmin"])):which.min(abs(scanseq-peakBounds[i,"scanmax"]))
                                   sapply(eicsRaw[[i]], function(eic) {
                                       f <- (2^eic(scanseq))-1
                                       return(sum(f[wh])*scanstep)
                                   })
                               }))
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) {
        message(sprintf("[getEICsAndQuantify]  .. done in %.1f secs.", stime))
    }
    return(quantmat)
}

slicepi <- function(object, cutoff = NULL, verbose = TRUE) {
    ## Set verbosity options
    subverbose <- max(as.integer(verbose) - 1L, 0)

    metadata <- list()

    ## If a cutoff is not supplied, compute it
    if (is.null(cutoff)) {
        if(verbose) {
            message("[slicepi] Computing cutoff")
        }
        cutoff <- getCutoff(object = object, mzSpacing = 2, verbose = subverbose)
        ## Get density quantiles and choose the higher of the two
        qcutoff <- which.min(abs(cutoff-densityQuantiles(object)))
        qs <- seq(0.001,0.999,0.001)
        qref <- which.min(abs(qs-0.99))
        metadata[["densityCutoff"]] <- max(cutoff, densityQuantiles(object)[qref])
    } else {
        metadata[["densityCutoff"]] <- cutoff
    }
    metadata[["densityQuantiles"]] <- densityQuantiles(object)
    
    if(verbose) {
        message("[slicepi] Computing peak bounds")
    }
    peakBounds <- computePeakBounds(densmat = densityEstimate(object), dcutoff = metadata[["densityCutoff"]], verbose = subverbose)

    ## Get EICs and quantifications
    if(verbose) {
        message("[slicepi] Quantifying peaks")
    }
    peakQuants <- getEICsAndQuantify(object = object, peakBounds = peakBounds, verbose = subverbose)

    ## Create SummarizedExperiment container
    CMSslice(assays = SimpleList(peakQuants = peakQuants),
             rowData = DataFrame(peakBounds),
             colData = colData(object),
             metadata = metadata,
             mzParams = .mzParams(object))
}
