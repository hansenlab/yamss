\name{getTIC}
\alias{getTIC}
\title{Compute total ion chromatogram (TIC) for a sample.}
\usage{
getTIC(object, sample)
}
\arguments{
\item{object}{An object of class \code{CMSraw}.}

\item{sample}{An integer - for which sample should the TIC be computed?.}
}
\value{
A vector with length equal to the number of scans containing the log2
sum of intensities at each scan.
}
\description{
Computes total ion chromatogram (TIC) for a single sample. Intensities 
are on the log2 scale. This requires a \code{CMSraw} object, typically 
produced from \code{readMSdata}.
}
\examples{
data(cmsRawExample)
tic <- getTIC(cmsRawExample, sample = 1)
head(tic)
}