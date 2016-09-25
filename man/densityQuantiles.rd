\name{densityQuantiles}
\alias{densityQuantiles}
\title{Get quantiles of the density estimate.}
\usage{
densityQuantiles(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A numeric vector containing the quantiles of the nonzero density values.
}
\description{
Gets the quantiles of the nonzero density values.
}

\examples{
cmsobj <- new("cms")
densityQuantiles(cmsobj)
}