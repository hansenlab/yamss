\name{slicepi}
\alias{slicepi}
\title{Process raw data to compute density estimate.}
\usage{
slicepi(cmsProc, cutoff = NULL, verbose = TRUE)
}
\arguments{
\item{cmsProc}{An object of class \code{CMSproc}.}

\item{cutoff}{A number indicating the threshold to apply to the density
estimate. \code{NULL} indicates that a data-driven threshold should be chosen.}

\item{verbose}{Should the function be verbose?}
}
\value{
An object of class \code{CMSslice} containing peak bounds and quantifications 
as well as sample and preprocessing metadata.
}
\description{
The \code{slicepi} method uses the bivariate approximate kernel density 
estimate computed by \code{bakedpi} and uses a cutoff to bound and quantify 
peaks.
}
\details{
\code{slicepi} uses the most intense features in set regions of the M/Z space 
to identify a data-driven density cutoff to detect peaks. Once peak bounds have 
been computed, the extracted ion chromatograms for the peaks are computed, and 
the EICs are integrated to obtain peak quantifications.
}
\examples{
densmat <- matrix(rnorm(600), nrow = 20, ncol = 30)
colnames(densmat) <- 1:ncol(densmat)
rownames(densmat) <- seq(350, by = 0.005, length.out = nrow(densmat))
densityQuantiles <- quantile(densmat, seq(from = 0, to = 1, by = 0.001))
cmsobj <- new("CMSproc", density = densmat, densityQuantiles = densityQuantiles)
cmsSlice <- slicepi(cmsobj, cutoff = densityQuantiles[999])
}