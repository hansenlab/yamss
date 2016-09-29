setOldClass("data.table")
setClass("CMSraw",
         slots = c(colData = "DataFrame",
                   rawDT = "data.table",
                   mzParams = "list"),
         prototype = prototype(
             rawDT = data.table(),
             colData = DataFrame(),
             mzParams = list()
         )
         )

setClass("CMSproc",
         contains = "CMSraw",
         representation( # FIXME: we might potentially want to store arguments of the call that made the object
             rtAlign = "logical", 
             bgcorrDT = "data.table",
             density = "matrix",
             densityQuantiles = "numeric"),
         prototype = prototype(
             rawDT = data.table(),
             colData = DataFrame(),
             mzParams = list(),
             rtAlign = FALSE,
             bgcorrDT = data.table(),
             density = matrix(),
             densityQuantiles = numeric()
         )
         )

setClass("CMSslice",
         contains = "SummarizedExperiment",
         representation(mzParams = "list")
         )

CMSslice <- function(..., mzParams = list()) {
    se <- SummarizedExperiment(...)
    cms <- as(se, "CMSslice")
    cms@mzParams <- mzParams
    cms
}

.printMZparams <- function(mzParams) {
    cat(sprintf("Number of scans: %i\n", mzParams$maxScan))
    cat(sprintf("M/Z: %f - %f\n", mzParams$minMZraw, mzParams$maxMZraw))
}
    
setMethod("show", signature(object = "CMSraw"),
          function(object) {
    cat(sprintf("class: %s\n", class(object)))
    cat(sprintf("Representing %i data files\n", nrow(object@colData)))
    .printMZparams(object@mzParams)
})

setMethod("show", signature(object = "CMSproc"),
          function(object) {
    cat(sprintf("class: %s\n", class(object)))
    cat(sprintf("Representing %i data files\n", nrow(object@colData)))
    .printMZparams(object@mzParams)
})



## Convenience functions
.minMZ <- function(object) {
    object@mzParams[["minMZ"]]
}

.maxMZ <- function(object) {
    object@mzParams[["maxMZ"]]
}

.maxScan <- function(object) {
    object@mzParams[["maxScan"]]
}

.sampleNumber <- function(object) {
    object@colData[, "sample"]
}

.rawDT <- function(object) {
    stopifnot(is(object, "CMSraw"))
    object@rawDT
}


## Exported accessors
setMethod("colData", signature(x = "CMSraw"), function(x) {
    x@colData
})

densityEstimate <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@density
}

densityCutoff <- function(object) {
    stopifnot(is(object, "CMSslice"))
    metadata(object)[["densityCutoff"]]
}

densityQuantiles <- function(object) {
    stopifnot(is(object, "CMSproc") | is(object, "CMSslice"))
    if (is(object, "CMSproc")) {
        return(object@densityQuantiles)
    }
    if (is(object, "CMSslice")) {
        return(metadata(object)[["densityQuantiles"]])
    }
}

peakBounds <- function(object) {
    stopifnot(is(object, "CMSslice"))
    rowData(object)
}

peakQuants <- function(object) {
    stopifnot(is(object, "CMSslice"))
    assay(object, "peakQuants")
}
