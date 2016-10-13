.getPeakInfo <- function(files) {
    isCDF <- length(grep("\\.cdf", files[1], ignore.case = TRUE))==1
    if (isCDF) {
        backend <- "netCDF"
    } else {
        backend <- "Ramp"
    }
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

.setMZParams <- function(rawDT) {
    mzParams <- list(
        minScan = min(rawDT[,scan]),
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
    rawDT <- .rawDT(object)
    setkey(rawDT, mz, scan)
    mzseq <- seq(as.integer(mzsubset[1]*1e5), as.integer(mzsubset[2]*1e5))
    .rawDT(object) <- rawDT[.(mzseq), nomatch = 0]
    .mzParams(object) <- .setMZParams(.rawDT(object))
    object
}

readMSdata <- function(files, colData = NULL,
                       mzsubset = NULL, verbose = FALSE) {
    if(verbose) {
        message(sprintf("[readRaw]: Reading %i files", length(files)))
    }
    stopifnot(all(file.exists(files)))
    if (!is.null(colData)) {
        stopifnot(is(colData, "DataFrame"))
        stopifnot(nrow(colData)==length(files))
        stopifnot(!any(c("sample", "files") %in% colnames(colData)))
    }
    cmsRaw <- new("CMSraw")
    rawPeakInfo <- .getPeakInfo(files)
    ## Make raw data matrix and data.table
    rawdatamat <- do.call(rbind,
                          lapply(seq_along(rawPeakInfo), function(s) {
                              cbind(do.call(rbind,
                                            lapply(seq_along(rawPeakInfo[[s]]), function(scan) {
                                                cbind(rawPeakInfo[[s]][[scan]], scan)
                                            })), s)
                          }))
    colnames(rawdatamat) <- c("mz", "intensity", "scan", "sample")
    rawDT <- data.table(mz = as.integer(rawdatamat[,"mz"]*1e5),
                     intensity = rawdatamat[,"intensity"],
                     scan = rawdatamat[,"scan"],
                     sample = rawdatamat[,"sample"])
    .mzParams(cmsRaw) <- .setMZParams(rawDT)
    .rawDT(cmsRaw) <- rawDT
    fileData <- DataFrame(sample = seq_along(files), files = files)
    if(is.null(colData)) {
        colData(cmsRaw) <- fileData
    } else {
        colData(cmsRaw) <- cbind(fileData, colData)
    }
    cmsRaw <- .subsetByMZ(cmsRaw, mzsubset = mzsubset)
    return(cmsRaw)
}
