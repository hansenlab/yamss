\name{files}
\alias{files}
\title{Get filepaths.}
\usage{
files(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A character vector of filepaths to the raw data.
}
\description{
Gets the raw data filepaths.
}

\examples{
cmsobj <- new("cms")
files(cmsobj)
}