\name{getTIC}
\alias{getTIC}
\title{Compute total ion chromatogram (TIC).}
\usage{
getTIC(object)
}
\arguments{
\item{object}{An object of class \code{CMSraw}.}
}
\value{
A vector with length equal to the number of scans containing the log2
sum of intensities at each scan.
}
\description{
Computes TIC from parsed raw data. Intensities are on the log2 scale. This
can be computed after running \code{getPeakInfo}.
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1", recursive = TRUE, full.names = TRUE)[1:2]
colData <- DataFrame(files = files)
cmsRaw <- readMSdata(files = files, colData = colData, verbose = TRUE)
tic <- getTIC(cmsRaw)
head(tic)
}
}