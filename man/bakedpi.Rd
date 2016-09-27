\name{bakedpi}
\alias{bakedpi}
\alias{updatePeaks}
\title{Process raw data to obtain quantifications.}
\usage{
bakedpi(files, dbandwidth = c(0.005, 10), dgridstep = c(0.005, 1),
  outfileDens = NULL, dortalign = FALSE, mzsubset = NULL, verbose = TRUE)
updatePeaks(cms, cutoff)
}
\arguments{
\item{files}{A character vector of filepaths to the raw data}

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

\item{cms}{An object of class \code{CMS}.}

\item{cutoff}{A number indicating the new cutoff to use to threshold 
the density estimate.}
}
\value{
An object of class \code{CMS} containing preprocessing information.
}
\description{
The \code{bakedpi} Performs background correction, retention time
correction, peak detection, and differential analysis. As part of this
analysis, a two-dimenstional density is computed and thresholded (at a
data-driven threshold).

The \code{updatePeaks} function allows a user to change the density
threshold as described in the previous paragraph, to a user-defined
threshold. After applying the threshold, peaks are detected and
quantified again.
}
\details{
\code{bakedpi} first performs region-specific background correction. An
optional retention time correction step follows in which M/Z region-specific
shifts are computed to align the raw data. Next the two-dimensional density
estimate is computed. A cutoff is determined to call peaks which are then
quantified and run through a differential analysis pipeline.
}
\examples{
\dontrun{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1", recursive = TRUE, full.names = TRUE)
cmsobj <- bakedpi(files = files, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                  outfileDens = NULL, dortalign = TRUE, mzsubset = c(500, 520))
}
}
}