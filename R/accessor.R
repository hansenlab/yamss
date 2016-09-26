setClass("CMS",
         slots = c(fileNames = "character",
                   rtalign = "logical",
                   rawPeakInfo = "list",
                   rawPeakDT = "data.frame",
                   bgcorrDT = "data.frame",
                   mzParams = "list",
                   bgSmooths = "list",
                   density = "matrix",
                   densityCutoff = "numeric",
                   densityQuantiles = "numeric",
                   alignments = "list",
                   eicsRaw = "list",
                   eicsImputed = "list",
                   peakBounds = "matrix",
                   peakQuants = "matrix")
         )

setMethod("show", signature(object = "CMS"),
          function(object) {
})

getFileNames <- function(obj) {
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
