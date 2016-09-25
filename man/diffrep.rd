\name{diffrep}
\alias{diffrep}
\title{Get differential analysis report.}
\usage{
diffrep(obj)
}
\arguments{
\item{obj}{An object of class "cms".}
}
\value{
A \code{data.frame} containing the results of a differential analysis.
}
\description{
Gets the differential analysis results table.
}

\examples{
cmsobj <- new("cms")
diffrep(cmsobj)
}