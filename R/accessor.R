setClass("CMSraw",
         slots = c(phenoData = "DataFrame",
                   rawPeakDT = "data.frame",
                   mzParams = "list"
                   )
         )

setClass("CMSproc",
         slots = c(raw = "CMSraw",
                   rtalign = "logical",
                   bgcorrDT = "data.frame",
                   density = "matrix",
                   densityCutoff = "numeric",
                   densityQuantiles = "numeric",
                   peakBounds = "matrix",
                   peakQuants = "matrix")
         )

setMethod("show", signature(object = "CMSraw"),
          function(object) {
            cat(sprintf("An object of class CMSraw containing %i files", length(obj@fileNames)))
            cat(sprintf("M/Z: %f - %f. %i scans.", object@mzParams$minMZraw, object@mzParams$maxMZraw, object@mzParams$maxScan))
})

fileNames <- function(obj) {
    obj@fileNames
}

densityEstimate <- function(obj) {
    obj@density
}

densityCutoff <- function(obj) {
    obj@densityCutoff
}

densityQuantiles <- function(obj) {
    obj@densityQuantiles
}

peakBounds <- function(obj) {
    obj@peakBounds
}

peakQuants <- function(obj) {
    obj@peakQuants
}
