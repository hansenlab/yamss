utils::globalVariables(c("peaknum", "mzmin", "mzmax", "weight", "bg", "gmz",
                         "mz", "intensity", "eic", "tic", "scanmax", 
                         "scanmin", "scanorig", "N", "."))


.isArgumentTwoVector <- function(vec) {
    stopifnot(is.numeric(vec))
    stopifnot(length(vec) == 2)
    stopifnot(vec[1] <= vec[2])
    stopifnot(vec[1] > 0)
}

.digestDataTableRaw <- function(dt, digits = 6) {
    content <- matrix("", nrow = nrow(dt), ncol = 4)
    colnames(content) <- c("mz", "intensity", "scan", "sample")
    content[, "mz"] <- as.character(dt$mz)
    content[, "scan"] <- as.character(dt$scan)
    content[, "sample"] <- as.character(dt$sample)
    content[, "intensity"] <- sprintf(paste0("%.", digits, "f"), dt$intensity)
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", digits)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}

.digestDataTableBG <- function(dt, digits = 6) {
    if(all(c("gmz", "gscan") %in% colnames(dt))) {
        content <- matrix("", nrow = nrow(dt), ncol = 6)
        colnames(content) <- c("mz", "intensity", "scan", "sample",
                               "gmz", "gscan")
        content[, "gmz"] <- sprintf(paste0("%.", digits, "f"), dt$gmz)
        content[, "gscan"] <- sprintf(paste0("%.", digits, "f"), dt$gscan)
    } else {
        content <- matrix("", nrow = nrow(dt), ncol = 4)
        colnames(content) <- c("mz", "intensity", "scan", "sample")
    }
    content[, "mz"] <- as.character(dt$mz)
    content[, "scan"] <- as.character(dt$scan)
    content[, "sample"] <- as.character(dt$sample)
    content[, "intensity"] <- sprintf(paste0("%.", digits, "f"), dt$intensity)
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", digits)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}
