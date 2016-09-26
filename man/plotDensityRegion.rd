\name{plotDensityRegion}
\alias{plotDensityRegion}
\title{Image plot of region of density estimate.}
\usage{
plotDensityRegion(obj, mzrange, scanrange)
}
\arguments{
\item{obj}{An object of class "cms".}

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
densmat <- matrix(rnorm(600), nrow = 20, ncol = 30)
colnames(densmat) <- 1:ncol(densmat)
rownames(densmat) <- seq(350, by = 0.005, length.out = nrow(densmat))
cmsobj <- new("cms", dens = densmat)
plotDensityRegion(cmsobj, mzrange = c(350.01, 350.03))
}