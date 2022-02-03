\name{readMSdata}
\alias{readMSdata}
\title{Read raw mass spectral data files.}
\usage{
readMSdata(files, colData, mzsubset, verbose)
}
\arguments{
\item{files}{A character vector of filenames pointing to the raw data.}

\item{colData}{A \code{DataFrame} of phenotype and sample information.}

\item{mzsubset}{A length-2 vector indicating a subset of the M/Z range to 
process. \code{NULL} otherwise.}

\item{verbose}{Should the function be verbose?}
}
\value{
A vector with length equal to the number of scans containing the log2
sum of intensities at each scan.
}
\description{
Creates a \code{CMSraw} object that contains a \code{data.table} of 
raw mass spectral information for all samples. The resulting object 
also stores phenotype and sample information. This object is the 
basic encapsulation of essential raw experimental data and serves as 
the output for further processing methods.
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzML")
file <- list.files(filepath, pattern = "MSpos-Ex1",
                   recursive = TRUE, full.names = TRUE)[1]
colData <- DataFrame(group = "wild-type")
cmsRaw <- readMSdata(files = file, colData = colData, verbose = TRUE)
}
}