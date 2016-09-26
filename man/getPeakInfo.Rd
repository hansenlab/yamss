\name{getPeakInfo}
\alias{getPeakInfo}
\title{Read raw data files.}
\usage{
getPeakInfo(files)
}
\arguments{
\item{files}{A character vector of raw data filenames.}
}
\value{
A list with length equal to the number of files containing M/Z, 
  scan, and intensity information. Each list element is itself a list of
  spectrum information.
}
\description{
Reads raw .mzdata and .cdf files and stores the data in a list.
}
\details{
Only MS1 information is kept, and records that have zero intensity are also
removed.
}
\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, recursive = TRUE, full.names = TRUE)
files <- files[c(1,2)]
pinfo <- getPeakInfo(files)
}
}