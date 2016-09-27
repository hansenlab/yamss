setOldClass("data.table")
setClass("CMSraw",
         slots = c(phenoData = "data.frame", # perhaps DataFrane, and probably a different name.
                   rawDT = "data.table",
                   mzParams = "list")
         )

setClass("CMSproc",
         contains = "CMSraw",
         representation( # FIXME: we might potentially want to store arguments of the call that made the object
             bgcorrDT = "data.table",
             density = "matrix",
             densityCutoff = "numeric",
             densityQuantiles = "numeric",
             peakBounds = "matrix",
             peakQuants = "matrix")
         )

setMethod("show", signature(object = "CMSraw"),
          function(object) {
    cat("An object of class 'CMSraw'\n")
    cat(sprintf("Representing %i data files\n", nrow(object@phenoData)))
    cat(sprintf("Number of scans: %i\n", object@mzParams$maxScan))
    cat(sprintf("M/Z: %f - %f\n", object@mzParams$minMZraw, object@mzParams$maxMZraw))
})

setMethod("show", signature(object = "CMSproc"),
          function(object) {
    cat("An object of class 'CMSproc'\n")
    cat(sprintf("Representing %i data files\n", nrow(object@phenoData)))
    cat(sprintf("Number of scans: %i\n", object@mzParams$maxScan))
    cat(sprintf("M/Z: %f - %f\n", object@mzParams$minMZraw, object@mzParams$maxMZraw))
    # FIXME some more information
})

## Convenience functions
.mzMin <- function(object) {
    object@mzParams[["mzMin"]]
}

.mzMax <- function(object) {
    object@mzParams[["mzMax"]]
}

.sampleNumber <- function(object) {
    object@phenoData[, "sample"]
}

## Accessors for CMSraw and CMSproc
phenoInfo <- function(object) {
    stopifnot(is(object, "CMSraw") | is(object, "CMSproc"))
    object@phenoData
}

rawDT <- function(object) {
    stopifnot(is(object, "CMSraw") | is(object, "CMSproc"))
    object@rawDT
}

## Accessors for CMSproc
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