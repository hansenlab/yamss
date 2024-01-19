utils::globalVariables(
    c(
        "peaknum", "mzmin", "mzmax", "weight", 
        "bg", "gmz", "mz", "intensity", "eic", "tic", 
        "scanmax", "scanmin", "scanorig", "N", "."
    )
)


.isArgumentTwoVector <- function(vec) {
    stopifnot(is.numeric(vec))
    stopifnot(length(vec) == 2)
    stopifnot(vec[1] <= vec[2])
    stopifnot(vec[1] > 0)
}

.digestDataTableRaw <- function(dt) {
    content <- matrix("", nrow = nrow(dt), ncol = 4)
    colnames(content) <- c("mz", "intensity", "scan", "sample")
    content[, "mz"] <- as.character(dt$mz)
    content[, "scan"] <- as.character(dt$scan)
    content[, "sample"] <- as.character(dt$sample)
    content[, "intensity"] <- sprintf("%.3f", dt$intensity)
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", 3)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}

.digestDataTableBG <- function(dt) {
    if (all(c("gmz", "gscan") %in% colnames(dt))) {
        content <- matrix("", nrow = nrow(dt), ncol = 6)
        colnames(content) <- c("mz", "intensity", "scan", "sample", "gmz", "gscan")
        content[, "gmz"] <- sprintf("%.3f", dt$gmz)
        content[, "gscan"] <- sprintf("%.3f", dt$gscan)
    } else {
        content <- matrix("", nrow = nrow(dt), ncol = 4)
        colnames(content) <- c("mz", "intensity", "scan", "sample")
    }
    content[, "mz"] <- as.character(dt$mz)
    content[, "scan"] <- as.character(dt$scan)
    content[, "sample"] <- as.character(dt$sample)
    content[, "intensity"] <- sprintf("%.3f", dt$intensity)
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", 3)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}

.digestPeakBounds <- function(rowData) {
    content <- matrix("", nrow = nrow(rowData), ncol = 5)
    colnames(content) <- c("mzmin", "mzmax", "scanmin", "scanmax", "peaknum")
    content[, "mzmin"] <- sprintf("%.2f", rowData[["mzmin"]])
    content[, "mzmax"] <- sprintf("%.2f", rowData[["mzmax"]])
    content[, "scanmin"] <- as.character(rowData[["scanmin"]])
    content[, "scanmax"] <- as.character(rowData[["scanmax"]])
    content[, "peaknum"] <- as.character(rowData[["peaknum"]])
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", 2)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}

.digestPeakQuants <- function(peakQuants) {
    content <- matrix(sprintf("%.1f", peakQuants), nrow = nrow(peakQuants), ncol = ncol(peakQuants))
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.0", rep("0", 1)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}
