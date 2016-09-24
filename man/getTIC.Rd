\name{getTIC}
\alias{getTIC}
\title{Compute total ion chromatogram (TIC).}
\usage{
getTIC(peak.info)
}
\arguments{
\item{peak.info}{A list of spectrum information for every scan.}
}
\value{
A vector with length equal to the number of scans (equal to the
  length of \code{peak.info}) containing the log2 sum of intensities
  at each scan.
}
\description{
Computes TIC from parsed raw data. Intensities are on the log2 scale. This
can be computed after running \code{getPeakInfo}.
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
file <- list.files(filepath, recursive = TRUE, full.names = TRUE)[1]
pinfo <- getPeakInfo(file)
tic <- getTIC(pinfo[[1]])
}
}

