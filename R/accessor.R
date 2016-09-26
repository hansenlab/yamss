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
                   densityQuants = "numeric",
                   alignments = "list",
                   xicsRaw = "list",
                   xicsImputed = "list",
                   peakBounds = "matrix",
                   peakQuants = "matrix")
         )

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
    obj@densityQuants
}

peakBounds <- function(obj) {
    obj@peakBounds
}

peakQuants <- function(obj) {
    obj@peakQuants
}
