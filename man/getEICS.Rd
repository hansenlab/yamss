\name{getEICS}
\alias{getEICS}
\title{Compute extracted ion chromatograms for multiple regions.}
\usage{
getEICS(cms, mzranges)
}
\arguments{
\item{cms}{An object of class \code{CMS}.}
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
\dontrun{
if (require(mtbls2)) {
data(mtbls2)
## Same example as in example(bakedpi) and the vignette
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1", recursive = TRUE, full.names = TRUE)
cmsobj <- bakedpi(files = files, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                  outfileDens = NULL, dortalign = TRUE, mzsubset = c(500, 520))
mzranges <- rbind(c(500.01, 500.03), c(511.44, 511.46))
eicList <- getEICS(cmsobj, mzranges)
}
}
}


