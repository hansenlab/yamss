\name{rawdataDT}
\alias{rawdataDT}
\title{Get raw spectral data in a data.table.}
\usage{
rawdataDT(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A \code{data.table} of the raw M/Z, scan, intensity, and sample 
information.
}
\description{
Gets the raw spectral information in \code{data.table} form.
}

\examples{
cmsobj <- new("cms")
rawdataDT(cmsobj)
}