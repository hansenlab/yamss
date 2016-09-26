\name{getDiffTable}
\alias{getDiffTable}
\title{Perform differential analysis.}
\usage{
getDiffTable(quants, classes)
}
\arguments{
\item{quants}{A quantification matrix with rows corresponding to features
(peaks) and columns corresponding to samples.}

\item{classes}{A character vector of class labels for the samples.}
}
\value{
A \code{data.frame} containing differential analysis information
  including log fold changes and p-values.
}
\description{
Performs differential abundance analysis on a peaks-by-samples matrix of
unlogged quantifications.
}
\details{
Differential analysis is performed using the \code{limma} package which
uses empirical Bayes methods in the estimation of feature-wise variances.
}
\examples{
quantmat <- matrix(rnorm(12*5000), nrow = 5000, ncol = 12)
classes <- rep(c("case", "control"), each = 6)
difftab <- getDiffTable(quantmat, classes)
}