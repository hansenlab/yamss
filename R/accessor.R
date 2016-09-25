files <- function(obj) {
    obj@files
}

classes <- function(obj) {
    obj@classes
}

rawpeakinfo <- function(obj) {
    obj@rawpeakinfo
}

densityEstimate <- function(obj) {
    obj@dens
}

densityCutoff <- function(obj) {
    obj@dcutoff
}

densityQuantiles <- function(obj) {
    obj@densquants
}

peakBounds <- function(obj) {
    obj@blobs
}

peakQuants <- function(obj) {
    obj@quants
}

diffrep <- function(obj) {
    obj@diffrep
}