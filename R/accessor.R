setClass("CMSraw",
         slots = c(phenoData = "DataFrame", # perhaps data.frame, and probably a different name.
                   rawDT = "data.frame",
                   mzParams = "list"
                   )
         )

setClass("CMSproc",
         contains = "CMSraw",
         representation( # FIXME: we might potentially want to store arguments of the call that made the object
             bgcorrDT = "data.frame",
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

.mzMin <- function(object) {
    object@mzParams[["mzMin"]]
}

.mzMax <- function(object) {
    object@mzParams[["mzMax"]]
}

.sampleNumber <- function(object) {
    object@phenoData[, "sample"]
}

## FIMXE: need some phenoData extractor


densityEstimate <- function(object) {
    object@density
}

densityCutoff <- function(object) {
    object@densityCutoff
}

densityQuantiles <- function(object) {
    object@densityQuantiles
}

peakBounds <- function(object) {
    object@peakBounds
}

peakQuants <- function(object) {
    object@peakQuants
}
