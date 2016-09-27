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
             densityCutoff = "numeric",
             densityQuantiles = "numeric",
             peakBounds = "matrix",
             peakQuants = "matrix"),
         prototype = prototype(
             rawDT = data.table(),
             colData = DataFrame(),
             mzParams = list(),
             rtAlign = FALSE,
             bgcorrDT = data.table(),
             density = matrix(),
             densityCutoff = numeric(),
             densityQuantiles = numeric(),
             peakBounds = matrix(),
             peakQuants = matrix()
         )
         )

setMethod("show", signature(object = "CMSraw"),
          function(object) {
    cat("An object of class 'CMSraw'\n")
    cat(sprintf("Representing %i data files\n", nrow(object@colData)))
    cat(sprintf("Number of scans: %i\n", object@mzParams$maxScan))
    cat(sprintf("M/Z: %f - %f\n", object@mzParams$minMZraw, object@mzParams$maxMZraw))
})

setMethod("show", signature(object = "CMSproc"),
          function(object) {
    cat("An object of class 'CMSproc'\n")
    cat(sprintf("Representing %i data files\n", nrow(object@colData)))
    cat(sprintf("Number of scans: %i\n", object@mzParams$maxScan))
    cat(sprintf("M/Z: %f - %f\n", object@mzParams$minMZraw, object@mzParams$maxMZraw))
    cat(sprintf("Number of peaks: %i\n", nrow(object@peakBounds)))
})

## Convenience functions
.mzMin <- function(object) {
    object@mzParams[["mzMin"]]
}

.mzMax <- function(object) {
    object@mzParams[["mzMax"]]
}

.sampleNumber <- function(object) {
    object@colData[, "sample"]
}

.rawDT <- function(object) {
    stopifnot(is(object, "CMSraw"))
    object@rawDT
}


## Accessors for CMSproc
setMethod("colData", signature(x = "CMSraw"), function(x) {
    x@colData
})

densityEstimate <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@density
}

densityCutoff <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@densityCutoff
}

densityQuantiles <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@densityQuantiles
}

peakBounds <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@peakBounds
}

peakQuants <- function(object) {
    stopifnot(is(object, "CMSproc"))
    object@peakQuants
}
