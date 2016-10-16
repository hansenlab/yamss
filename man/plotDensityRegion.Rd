\name{plotDensityRegion}
\alias{plotDensityRegion}
\title{Image plot of region of density estimate.}
\usage{
plotDensityRegion(cms, mzrange, scanrange)
}
\arguments{
\item{cms}{An object of class \code{CMSproc}.}

\item{mzrange}{A length-2 vector indicating the M/Z range to plot.}

\item{scanrange}{A length-2 vector indicating the scan range to plot.}
}
\description{
Makes an \code{image} plot of the density estimate in the specified 
M/Z and scan region.
}
\value{
  This function is invoked for its side effect of plotting.
}
\examples{
## For illustration purposes, we make a "dummy" object
## with a random matrix as the density estimate

densmat <- matrix(rnorm(600), nrow = 20, ncol = 30)
colnames(densmat) <- 1:ncol(densmat)
rownames(densmat) <- seq(350, by = 0.005, length.out = nrow(densmat))
densityQuantiles <- quantile(densmat, seq(from = 0, to = 1, by = 0.001))
cmsobj <- new("CMSproc", density = densmat, densityQuantiles = densityQuantiles)
plotDensityRegion(cmsobj, mzrange = c(350.01, 350.03), scanrange = c(10,20))
}