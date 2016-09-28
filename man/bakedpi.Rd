\name{bakedpi}
\alias{bakedpi}
\title{Process raw data to compute density estimate.}
\usage{
bakedpi(cmsRaw, dbandwidth = c(0.005, 10), dgridstep = c(0.005, 1),
  outfileDens = NULL, dortalign = FALSE, mzsubset = NULL, verbose = TRUE)
}
\arguments{
\item{cmsRaw}{An object of class \code{CMSraw}.}

\item{dbandwidth}{A length-2 vector indicating the kernel density bandwidth
in the M/Z and retention time (scan) directions. Default: \code{c(0.005,10)}}

\item{dgridstep}{A length-2 vector indicating the grid step sizes. Default: 
\code{c(0.005,1)}.}

\item{outfileDens}{Name of a file to save density estimate. If NULL,
no output is saved.}

\item{dortalign}{A logical value. Should retention time correction be
performed?}

\item{mzsubset}{A length-2 vector indicating a subset of the M/Z range to 
process. \code{NULL} otherwise.}

\item{verbose}{Should the function be verbose?}
}
\value{
An object of class \code{CMSproc} containing background corrected intensities, 
the bivariate kernel density estimate, and quantiles of the nonzero values in 
the density estimate.
}
\description{
The \code{bakedpi} method stands for bivariate approximate kernel density 
estimation for peak identification. It performs background correction, 
retention time correction, and bivariate kernel density estimation.
}
\details{
\code{bakedpi} first performs region-specific background correction. An
optional retention time correction step follows in which M/Z region-specific
shifts are computed to align the raw data. Next the two-dimensional density
estimate is computed. The purpose of this function is to take the raw data 
read in by \code{readMSdata} and perform the steps necessary for bivariate 
kernel density estimation. The output of this function is used by 
\code{slicepi} to detect peaks and provide peak quantifications.
}
\examples{
\dontrun{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1", recursive = TRUE, full.names = TRUE)
colData <- DataFrame(sampClasses = rep(c("wild-type", "mutant"), each = 4))
cmsRaw <- readMSdata(files = files, colData = colData, verbose = TRUE)
cmsProc <- bakedpi(cmsRaw, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                   outfileDens = NULL, dortalign = TRUE, mzsubset = c(500, 520))
}
}
}