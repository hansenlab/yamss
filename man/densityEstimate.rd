\name{densityEstimate}
\alias{densityEstimate}
\title{Get kernel density estimate.}
\usage{
densityEstimate(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A matrix containing the kernel density estimate with rows 
corresponding to M/Z values and columns corresponding to scans.
}
\description{
Gets the density estimate matrix.
}

\examples{
cmsobj <- new("cms")
densityEstimate(cmsobj)
}