utils::globalVariables(c("blobnum", "mzmin", "mzmax", "weight", "bg", "gmz",
                         "mz", "intensity",
                         "xic", "scanmax", "scanmin", "scanorig", "N", "."))


.isArgumentTwoVector <- function(vec) {
    stopifnot(is.numeric(vec))
    stopifnot(length(vec) == 2)
    stopifnot(vec[1] <= vec[2])
    stopifnot(vec[1] > 0)
}
