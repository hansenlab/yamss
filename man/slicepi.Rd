\name{slicepi}
\alias{slicepi}
\title{Process raw data to compute density estimate.}
\usage{
slicepi(object, cutoff = NULL, verbose = TRUE)
}
\arguments{
\item{object}{An object of class \code{CMSproc}.}

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
data(cmsRawExample)
cmsProc <- bakedpi(cmsRawExample, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                   outfileDens = NULL, dortalign = FALSE, verbose = TRUE)
dqs <- densityQuantiles(cmsProc)
cmsSlice <- slicepi(cmsProc, cutoff = dqs[996], verbose = TRUE)
cmsSlice
}
