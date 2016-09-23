#' Read raw data files.
#' 
#' Reads raw .mzdata and .cdf files and stores the data in a \code{data.table}.
#' 
#' Only MS1 information is kept, and records that have zero intensity are also
#' removed.
#' 
#' @param files A character vector of raw data filenames.
#' @return A list with length equal to the number of files containing M/Z, 
#'   scan, and intensity information. Each list element is itself a list of
#'   spectrum information.
#' 
#' @examples
#' \dontrun{
#' if (require(mtbls2)) {
#' data(mtbls2)
#' filepath <- file.path(find.package("mtbls2"), "mzData")
#' files <- list.files(filepath, recursive = TRUE, full.names = TRUE)
#' files <- files[c(1,2)]
#' pinfo <- getPeakInfo(files)
#' }
#' }
getPeakInfo <- function(files) {
    peakInfoAllfiles <- lapply(seq_along(files), function(i) {
        msobj <- openMSfile(files[i])
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

#' Compute total ion chromatogram (TIC).
#' 
#' Computes TIC from parsed raw data. Intensities are on the log2 scale. This
#' can be computed after running \code{getPeakInfo}.
#' 
#' @param peak.info A list of spectrum information for every scan.
#' @return A vector with length equal to the number of scans (equal to the
#'   length of \code{peak.info}) containing the log2 sum of intensities
#'   at each scan.
#' 
#' @examples
#' \dontrun{
#' if (require(mtbls2)) {
#' data(mtbls2)
#' filepath <- file.path(find.package("mtbls2"), "mzData")
#' file <- list.files(filepath, recursive = TRUE, full.names = TRUE)[1]
#' pinfo <- getPeakInfo(file)
#' tic <- getTIC(pinfo[[1]])
#' }
#' }
getTIC <- function(peak.info) {
    log2(sapply(peak.info, function(spec) { sum(spec[,2]) }) + 1)
}

#' Compute extracted ion chromatograms for multiple regions.
#' 
#' Computes XICs for the given M/Z ranges. Intensities are on the log2 scale.
#' 
#' In a given M/Z range, the maximum intensity observed in each scan gives the
#' extracted ion chromatogram.
#' 
#' @param mzranges A 2-column matrix where each row corresponds to one M/Z
#'   range and the first and second columns are the minimum and maximum
#'   M/Z values for the range respectively.
#' @param DT A \code{data.table} (usually from \code{readRawDataAsDataTable})
#'   containing M/Z, scan, intensity, and sample information.
#' @return A list with length equal to the number of rows of \code{mzranges}
#'   where each list element is a # scans by # samples matrix of XICs (on the
#'   log2 scale).
#' 
#' @examples
#' \dontrun{
#' if (require(mtbls2) & require(yamss)) {
#' data(mtbls2)
#' filepath <- file.path(find.package("mtbls2"), "mzData")
#' files <- list.files(filepath, recursive = TRUE, full.names = TRUE)
#' obj <- new("cms")
#' obj@files <- files
#' ## Parse raw data
#' out <- readRawDataAsDataTable(obj = obj, files = files)
#' obj <- out$obj
#' DT <- out$DT
#' mzranges <- rbind(c(455.09, 455.13), c(264, 264.04))
#' xicList <- getAllXics(mzranges, DT)
#' }
#' }
getAllXics <- function(mzranges, DT) {
    ## Convert mzranges to an IRanges
    mzranges <- IRanges(start = as.integer(mzranges[,1]*1e5), end = as.integer(mzranges[,2]*1e5))
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

#' Perform differential analysis.
#' 
#' Performs differential abundance analysis on a peaks-by-samples matrix of
#' unlogged quantifications.
#' 
#' Differential analysis is performed using the \code{limma} package which
#' uses empirical Bayes methods in the estimation of feature-wise variances.
#' 
#' @param quants A quantification matrix with rows corresponding to features
#'   (peaks) and columns corresponding to samples.
#' @param classes A character vector of class labels for the samples.
#' @return A \code{data.frame} containing differential analysis information
#'   including log fold changes and p-values.
#' 
#' @examples
#' quantmat <- matrix(rnorm(12*5000), nrow = 5000, ncol = 12)
#' classes <- rep(c("case", "control"), each = 6)
#' difftab <- getDiffTable(quantmat, classes)
getDiffTable <- function(quants, classes) {
    df <- data.frame(classes = classes)
    design <- model.matrix(~classes, data = df)
    fit <- lmFit(log2(quants + 1), design = design)
    fit <- eBayes(fit)
    return(topTable(fit, coef = 2, number = Inf, sort.by = "P"))
}
