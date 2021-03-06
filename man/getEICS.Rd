\name{getEICS}
\alias{getEICS}
\title{Compute extracted ion chromatograms for multiple regions.}
\usage{
getEICS(object, mzranges)
}
\arguments{
\item{object}{An object of class \code{CMSraw} or \code{CMSproc}.}
\item{mzranges}{A 2-column matrix where each row corresponds to one M/Z
range and the first and second columns are the minimum and maximum
M/Z values for the range respectively.}
}
\value{
A list with length equal to the number of rows of \code{mzranges}
  where each list element is a # scans by # samples matrix of EICs (on the
  log2 scale).
}
\description{
Computes extracted ion chromatrograms (EICs) for the given M/Z ranges. 
Intensities are on the log2 scale.
}
\details{
In a given M/Z range, the maximum intensity observed in each scan gives the
extracted ion chromatogram.
}
\examples{
data(cmsRawExample)
mzranges <- rbind(c(500.01, 500.03), c(501.3, 501.5))
eicList <- getEICS(cmsRawExample, mzranges)
}


