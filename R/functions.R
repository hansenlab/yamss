getTIC <- function(object, sample) {
    stopifnot(is(object, "CMSraw"))
    rawDT <- .rawDT(object)
    setkey(rawDT, scan)
    ticDT <- rawDT[sample==sample, .(tic = log2(sum(intensity)+1)), by = scan]
    tic <- rep(0, .maxScan(object))
    tic[ticDT[,scan]] <- ticDT[,tic]
    return(tic)
}

getEICS <- function(object, mzranges) {
    stopifnot(is(object, "CMSraw") | is(object, "CMSproc"))
    ## Convert mzranges to an IRanges if a matrix
    if (class(mzranges) != "IRanges") {
        mzranges <- IRanges(start = as.integer(mzranges[,1]*1e5),
                            end = as.integer(mzranges[,2]*1e5))
    }
    ## Sort DT and get XICs
    rawDT <- object@rawDT
    setkey(rawDT, mz, scan, sample)
    maxScan <- object@mzParams$maxScan
    numSamples <- length(unique(rawDT[,sample]))
    eics <- lapply(1:length(mzranges), function(i) {
        x <- rawDT[.(seq(start(mzranges[i]), end(mzranges[i]))), nomatch = 0]
        if (nrow(x)==0)
            return(matrix(0, nrow = maxScan, ncol = numSamples))
        x <- x[, eic := log2(max(intensity)+1), by = .(scan, sample)]
        dup <- duplicated(x[,.(scan,sample)])
        x <- x[!dup]
        m <- as.matrix(sparseMatrix(i = x[,scan], j = x[,sample], x = x[,eic]))
        if (nrow(m) < maxScan) {
            m <- rbind(m, matrix(0, nrow = maxScan-nrow(m), ncol = ncol(m)))
        }
        if (ncol(m) < numSamples) {
            m <- cbind(m, matrix(0, nrow = nrow(m), ncol = numSamples-ncol(m)))
        }
        return(m)
    })
    return(eics)
}

diffrep <- function(cms, classes) {
    stopifnot(is(cms, "CMSproc"))
    quants <- cms@peakQuants
    df <- data.frame(classes = classes)
    design <- model.matrix(~classes, data = df)
    fit <- lmFit(log2(quants + 1), design = design)
    fit <- eBayes(fit)
    out <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
    out
}

plotDensityRegion <- function(cms, mzrange, scanrange) {
    stopifnot(is(cms, "CMSproc"))
    .isArgumentTwoVector(mzrange)
    .isArgumentTwoVector(scanrange)
    if (nrow(cms@density)==0) {
        stop("'CMS' cmsect must have a density estimate")
    }
    mzs <- as.numeric(rownames(cms@density))
    scans <- as.numeric(colnames(cms@density))
    idxMZ <- which.min(abs(mzrange[1]-mzs)):which.min(abs(mzrange[2]-mzs))
    idxScan <- which.min(abs(scanrange[1]-scans)):which.min(abs(scanrange[2]-scans))
    subdensmat <- cms@density[idxMZ, idxScan]

    mypalette <- colorRampPalette(c("white", "palegoldenrod",
                                    "palegreen", "#99ccff", "#ff9999", "red"))
    colorsdens <- c(rep("white", 890), mypalette(110))
    main <- sprintf("M/Z: %f - %f. Scan: %i - %i", mzrange[1], mzrange[2],
                    scanrange[1], scanrange[2])
    image(z = t(subdensmat), x = scanrange[1]:scanrange[2], y = mzs[idxMZ],
          col = colorsdens, breaks = cms@densityQuantiles,
          xlab = "Scan", ylab = "M/Z", main = main)
}

updatePeaks <- function(cms, cutoff) {
    stopifnot(is(cms, "CMSproc"))
    if (nrow(cms@density)==0) {
        stop("cms cmsect must have a density estimate")
    }
    cms@peakBounds <- computePeakBounds(cms@density, dcutoff = cutoff)
    if (cms@rtAlign) {
        cms <- getXICsAndQuantifyWithRetentionTime(object = cms)
    } else {
        cms <- getXICsAndQuantifyWithoutRetentionTime(object = cms)
    }
    return(cms)
}
