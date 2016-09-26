\name{updatePeaks}
\alias{updatePeaks}
\title{Re-detect peaks and obtain quantifications.}
\usage{
updatePeaks(obj, cutoff)
}
\arguments{
\item{obj}{An object of class "cms".}

\item{cutoff}{A number indicating the new cutoff to use to threshold 
the density estimate.}
}
\value{
A "cms" object with updated peak bounds and quantifications corresponding 
to \code{cutoff}.
}
\description{
Applies a new cutoff to the density estimate to call peaks and re-quantifies 
the peaks.
}
\examples{
\dontrun{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1" recursive = TRUE, full.names = TRUE)
classes <- rep(c("wild-type", "mutant"), each = 4)
cmsobj <- bakedpi(files = files, classes = classes, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1), outfileDens = NULL, dortalign = TRUE, mzsubset = c(500, 520))
cmsobj2 <- updatePeaks(cmsobj, densityCutoff(cmsobj)*0.99)
}
}
}