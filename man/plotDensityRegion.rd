\name{plotDensityRegion}
\alias{plotDensityRegion}
\title{Image plot of region of density estimate.}
\usage{
plotDensityRegion(cms, mzrange, scanrange)
}
\arguments{
\item{cms}{An object of class \code{CMS}.}

\item{mzrange}{A length-2 vector indicating the M/Z range to plot.}

\item{scanrange}{A length-2 vector indicating the scan range to plot.}
}
\description{
Makes an \code{image} plot of the density estimate in the specified 
M/Z and scan region.
}

\examples{
## For illustration purposes, we make a "dummy" object
## with a random matrix as the density estimate
# FIXME
#densmat <- matrix(rnorm(600), nrow = 20, ncol = 30)
#colnames(densmat) <- 1:ncol(densmat)
#rownames(densmat) <- seq(350, by = 0.005, length.out = nrow(densmat))
#cmsobj <- new("CMS", density = densmat)
#plotDensityRegion(cmsobj, mzrange = c(350.01, 350.03), scanrange = c(10,20))
}