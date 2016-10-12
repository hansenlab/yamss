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

setMethod("show", signature(object = "CMSslice"),
          function(object) {
    ## FIXME add callNextMethod
    cat("An object of class 'CMSslice'\n")
    cat(sprintf("Representing %i data files\n", nrow(object@colData)))
    cat(sprintf("Number of scans: %i\n", object@mzParams$maxScan))
    cat(sprintf("M/Z: %f - %f\n", object@mzParams$minMZraw, object@mzParams$maxMZraw))
    cat(sprintf("Number of peaks: %i\n", nrow(object)))
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

.mzParams <- function(object) {
    object@mzParams
}

".mzParams<-" <- function(object, value) {
    object@mzParams <- value
    object
}

.sampleNumber <- function(object) {
    object@colData[, "sample"]
}

.rawDT <- function(object) {
    stopifnot(is(object, "CMSraw"))
    object@rawDT
}

".rawDT<-" <- function(object, value) {
    stopifnot(is(object, "CMSraw"))
    object@rawDT <- value
    object
}

.bgcorrDT <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@bgcorrDT
}

".bgcorrDT<-" <- function(object, value) {
    stopifnot(is(object, "CMSproc"))
    object@bgcorrDT <- value
    object
}

".rtAlign<-" <- function(object, value) {
    stopifnot(is(object, "CMSproc"))
    object@rtAlign <- value
    object
}

".densityEstimate<-" <- function(object, value) {
    stopifnot(is(object, "CMSproc"))
    object@density <- value
    object
}

".densityQuantiles<-" <- function(object, value) {
    stopifnot(is(object, "CMSproc"))
    object@densityQuantiles <- value
    object
}

## Exported accessors
setMethod("colData", signature(x = "CMSraw"), function(x) {
    x@colData
})

setGeneric("colData<-", function(object, value) standardGeneric("colData<-"))

setReplaceMethod("colData", "CMSraw", function(object, value) {
    object@colData <- value
    object
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
