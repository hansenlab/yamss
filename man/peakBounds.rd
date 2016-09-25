\name{peakBounds}
\alias{peakBounds}
\title{Get M/Z and scan bounds of called peaks.}
\usage{
peakBounds(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A matrix of information about the peaks detected: M/Z bounds, 
scan bounds, and ID number ("blobnum").
}
\description{
Gets the information on M/Z and scan bounds for detected peaks.
}

\examples{
cmsobj <- new("cms")
peakBounds(cmsobj)
}