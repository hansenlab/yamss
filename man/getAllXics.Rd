\name{getAllXics}
\alias{getAllXics}
\title{Compute extracted ion chromatograms for multiple regions.}
\usage{
getAllXics(mzranges, DT)
}
\arguments{
\item{mzranges}{A 2-column matrix where each row corresponds to one M/Z
range and the first and second columns are the minimum and maximum
M/Z values for the range respectively.}

\item{DT}{A \code{data.table} (usually from \code{readRawDataAsDataTable})
containing M/Z, scan, intensity, and sample information.}
}
\value{
A list with length equal to the number of rows of \code{mzranges}
  where each list element is a # scans by # samples matrix of XICs (on the
  log2 scale).
}
\description{
Computes XICs for the given M/Z ranges. Intensities are on the log2 scale.
}
\details{
In a given M/Z range, the maximum intensity observed in each scan gives the
extracted ion chromatogram.
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, recursive = TRUE, full.names = TRUE)
obj <- new("cms", files = files)
## Parse raw data
out <- readRawDataAsDataTable(obj = obj, files = files)
obj <- out$obj
DT <- out$DT
mzranges <- rbind(c(455.09, 455.13), c(264, 264.04))
xicList <- getAllXics(mzranges, DT)
}
}

