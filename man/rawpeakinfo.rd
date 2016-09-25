\name{rawpeakinfo}
\alias{rawpeakinfo}
\title{Get raw spectral data.}
\usage{
rawpeakinfo(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A list of raw spectral information for each scan in each sample.
}
\description{
Gets the raw spectral information.
}

\examples{
cmsobj <- new("cms")
rawpeakinfo(cmsobj)
}