\name{bgcorrdataDT}
\alias{bgcorrdataDT}
\title{Get background-corrected data in a data.table.}
\usage{
bgcorrdataDT(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A \code{data.table} of the background-corrected M/Z, scan, 
intensity, and sample information.
}
\description{
Gets the background-corrected data in \code{data.table} form.
}

\examples{
cmsobj <- new("cms")
bgcorrdataDT(cmsobj)
}