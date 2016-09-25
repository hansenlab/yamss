\name{classes}
\alias{classes}
\title{Get sample classes.}
\usage{
classes(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A character vector of class labels for the samples.
}
\description{
Gets the sample class labels.
}

\examples{
cmsobj <- new("cms")
classes(cmsobj)
}