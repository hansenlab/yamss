\name{getTIC}
\alias{getTIC}
\title{Compute total ion chromatogram (TIC) for a sample.}
\usage{
getTIC(object, sample)
}
\arguments{
\item{object}{An object of class \code{CMSraw}.}

\item{sample}{An integer - for which sample should the TIC be computed?.}
}
\value{
A vector with length equal to the number of scans containing the log2
sum of intensities at each scan.
}
\description{
Computes total ion chromatogram (TIC) for a single sample. Intensities 
are on the log2 scale. This requires a .
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
file <- list.files(filepath, pattern = "MSpos-Ex1", recursive = TRUE, full.names = TRUE)[1]
colData <- DataFrame(group = "wild-type")
cmsRaw <- readMSdata(files = file, colData = colData, verbose = TRUE)
tic <- getTIC(cmsRaw, sample = 1)
head(tic)
}
}