\name{readRawDataAsDataTable}
\alias{readRawDataAsDataTable}
\title{Compute extracted ion chromatograms for multiple regions.}
\usage{
readRawDataAsDataTable(obj, files, verbose = FALSE)
}
\arguments{
\item{obj}{An object of class "cms".}

\item{files}{A character vector of raw data filenames.}

\item{verbose}{Should the function be verbose?}
}
\value{
A list with elements \code{obj} and \code{DT} containing the updated "cms"
  object and \code{data.table} of M/Z, scan, intensity, and sample information.
}
\description{
Reads raw .mzdata and .cdf files and stores the data in both a list within the 
  updated "cms" object and also in a \code{data.table}.
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, recursive = TRUE, full.names = TRUE)[1:2]
obj <- new("cms", files = files)
## Parse raw data
out <- readRawDataAsDataTable(obj = obj, files = files)
obj <- out$obj
DT <- out$DT
}
}
