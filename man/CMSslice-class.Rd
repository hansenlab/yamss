\name{CMSslice-class}
\docType{class}
\alias{CMSslice-class}
\alias{CMSslice}
\alias{densityCutoff}
\alias{peakBounds}
\alias{peakQuants}
\alias{show,CMSslice-method}
\title{A class to hold chromatography-mass spectrometry peak bounds and
  quantifications.}

\description{
This class is based on the \code{SummarizedExperiment} class. It holds 
information on peak quantifications, M/Z and scan bounds, sample information, 
and preprocessing metadata. Objects of the class can be constructed
using \code{CMSslice}.
}

\section{Utility functions}{
  We have the following utility functions:

  \describe{
    \item{\code{show}:}{
      The show method; prints the object.
    }
  }
}

\section{Accessors}{
  We have the following accessor functions:

  \describe{
    \item{\code{colData}:}{
      Gets the \code{DataFrame} containing phenotype and sample information.
    }
    \item{\code{densityCutoff}:}{
      Gets the value used to threshold the density for peak calling.
    }
    \item{\code{densityQuantiles}:}{
      Gets the quantiles of the nonzero values in the density estimate.
    }
    \item{\code{peakBounds}:}{
      Gets the \code{DataFrame} of M/Z bounds, scan bounds, and ID numbers 
      for detected peaks.
    }
    \item{\code{peakQuants}:}{
      Gets the matrix of peak quantifications (rows: peaks, columns: samples).
    }
  }
}

\examples{
## Construct a fake class
quants <- matrix(rnorm(12*5000), nrow = 5000, ncol = 12)
bounds <- cbind(mzmin = seq(from = 100, to = 1100, length.out = 5000),
                mzmax = seq(from = 100, to = 1100, length.out = 5000) + 0.1,
                scan.min = rep(10,5000), scan.max = rep(20, 5000),
                peaknum = 1:5000)
cmsobj <- CMSslice(assays = SimpleList(peakQuants = quants),
                   rowData = DataFrame(bounds))
head(peakQuants(cmsobj))

## A better example which takes 20s to run
\dontrun{

data(cmsRawExample)
cmsProc <- bakedpi(cmsRawExample,
                   dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                   dortalign = TRUE, mzsubset = c(500,510))
cutoff <- densityQuantiles(cmsProc)["99.9%"]
sliced <- slicepi(cmsProc, cutoff = cutoff, verbose = TRUE)
sliced
}
}