getPeakInfo <- function(files) {
    isCDF <- length(grep("\\.cdf", files[1], ignore.case = TRUE))==1
    backend <- ifelse(isCDF, "netCDF", "Ramp")
    peakInfoAllfiles <- lapply(seq_along(files), function(i) {
        msobj <- openMSfile(files[i], backend = backend)
        peakInfo <- peaks(msobj)
        ## Only keep scans where the MZ level is 1 (MS versus MS-MS)
        headerInfo <- header(msobj)
        whMS1 <- which(headerInfo$msLevel==1)
        peakInfo <- peakInfo[whMS1]
                                        # Remove rows with zero intensity
        peakInfo <- lapply(peakInfo, function(spectrum) {
            keep <- spectrum[,2] > 1e-6
            return(spectrum[keep,,drop = FALSE])
        })
        close(msobj)
                                        # Store retention time information
        attr(peakInfo, "rt") <- headerInfo$retentionTime[whMS1]
        return(peakInfo)
    })
    return(peakInfoAllfiles)
}

getTIC <- function(peak.info) {
    log2(sapply(peak.info, function(spec) { sum(spec[,2]) }) + 1)
}

getAllXics <- function(mzranges, DT) {
    ## Convert mzranges to an IRanges if a matrix
    if (class(mzranges) != "IRanges") {
        mzranges <- IRanges(start = as.integer(mzranges[,1]*1e5), 
            end = as.integer(mzranges[,2]*1e5))
    }
    ## Sort DT and get XICs
    setkey(DT, mz, scan, sample)
    maxScan <- max(DT[,scan])
    numSamples <- length(unique(DT[,sample]))
    xics <- lapply(1:length(mzranges), function(i) {
        x <- DT[.(seq(start(mzranges[i]), end(mzranges[i]))), nomatch = 0]
        if (nrow(x)==0)
            return(matrix(0, nrow = maxScan, ncol = numSamples))
        x <- x[, xic := log2(max(intensity)+1), by = .(scan, sample)]
        dup <- duplicated(x[,.(scan,sample)])
        x <- x[!dup]
        m <- as.matrix(sparseMatrix(i = x[,scan], j = x[,sample], x = x[,xic]))
        if (nrow(m) < maxScan) {
            m <- rbind(m, matrix(0, nrow = maxScan-nrow(m), ncol = ncol(m)))
        }
        if (ncol(m) < numSamples) {
            m <- cbind(m, matrix(0, nrow = nrow(m), ncol = numSamples-ncol(m)))
        }
        return(m)
    })
    return(xics)
}

getDiffTable <- function(quants, classes) {
    df <- data.frame(classes = classes)
    design <- model.matrix(~classes, data = df)
    fit <- lmFit(log2(quants + 1), design = design)
    fit <- eBayes(fit)
    return(topTable(fit, coef = 2, number = Inf, sort.by = "P"))
}

plotDensityRegion <- function(obj, mzrange, scanrange) {
    mzs <- as.numeric(rownames(obj@dens))
    scans <- as.numeric(colnames(obj@dens))
    idxMZ <- which.min(abs(mzrange[1]-mzs)):which.min(abs(mzrange[2]-mzs))
    idxScan <- which.min(abs(scanrange[1]-scans)):which.min(abs(scanrange[2]-scans))
    subdensmat <- obj@dens[idxMZ, idxScan]

    mypalette <- colorRampPalette(c("white", "palegoldenrod", "palegreen", "#99ccff", "#ff9999", "red"))
    colorsdens <- c(rep("white", 890), mypalette(110))
    image(z = t(subdensmat), x = scanrange[1]:scanrange[2], y = mzs[idxMZ], col = colorsdens, breaks = obj@densquants, xlab = "Scan", ylab = "M/Z", main = paste0("M/Z: ", mzrange[1], " - ", mzrange[2], ". Scans: ", scanrange[1], " - ", scanrange[2]))
}

updatePeaks <- function(obj, cutoff) {
    newblobs <- getBlobs(obj@dens, dcutoff = cutoff, verbose = FALSE)
    obj@blobs <- newblobs
    if (obj@rtalign) {
        obj <- getXICsAndQuantifyWithRetentionTime(obj = obj, verbose = FALSE)
    } else {
        obj <- getXICsAndQuantifyWithoutRetentionTime(obj = obj, verbose = FALSE)
    }
    return(obj)
}