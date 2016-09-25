\name{peakQuants}
\alias{peakQuants}
\title{Get peak quantification matrix.}
\usage{
peakQuants(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A matrix of quantifications where rows correspond to peaks and 
columns correspond to samples.
}
\description{
Gets the peak quantification matrix.
}

\examples{
cmsobj <- new("cms")
peakQuants(cmsobj)
}