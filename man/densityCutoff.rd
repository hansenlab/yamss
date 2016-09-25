\name{densityCutoff}
\alias{densityCutoff}
\title{Get density cutoff.}
\usage{
densityCutoff(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A number indicating the density cutoff used to bound peak regions.
}
\description{
Gets the density cutoff value.
}

\examples{
cmsobj <- new("cms")
densityCutoff(cmsobj)
}