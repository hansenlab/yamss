\name{diffrep}
\alias{diffrep}
\title{Perform differential analysis.}
\usage{
diffrep(cms, classes)
}
\arguments{
\item{cms}{An object of class \code{CMS}.}

\item{classes}{A character vector of class labels for the samples.}
}
\value{
A \code{data.frame} containing differential analysis information
  including log fold changes and p-values.
}
\description{
Performs differential abundance analysis on quantification information 
in a \code{CMS} object.
}
\details{
Differential analysis is performed using the \code{limma} package which
uses empirical Bayes methods in the estimation of feature-wise variances.
}
\examples{
quantmat <- matrix(rnorm(12*5000), nrow = 5000, ncol = 12)
classes <- rep(c("case", "control"), each = 6)
cmsobj <- new("CMS", peakQuants = quantmat)
difftab <- diffrep(cmsobj, classes)
}