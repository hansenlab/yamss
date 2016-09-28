\name{CMSslice-class}
\docType{class}
\alias{CMSslice-class}
\alias{densityCutoff}
\alias{densityQuantiles}
\alias{peakBounds}
\alias{peakQuants}
\alias{show,CMSslice-method}
\title{A class to hold chromatography-mass spectrometry peak bounds and quantifications.}
\description{
This class is a wrapper for the \code{SummarizedExperiment} class. It holds 
information on peak quantifications, M/Z and scan bounds, sample information, 
and preprocessing metadata.
}

\section{Slots}{
  \describe{
    \item{\code{assays}:}{
      a list with a single element "peakQuants" that holds the feature by 
      samples matrix of peak quantifications.
    }
    \item{\code{rowData}:}{
      a \code{DataFrame} with rows corresponding to peaks containing peak 
      M/Z and scan (retention time) bound information.
    }
    \item{\code{colData}:}{
      a \code{DataFrame} of phenotype and sample information.
    }
    \item{\code{metadata}:}{
      a list with elements "densityCutoff" and "densityQuantiles" which hold 
      the threshold used to slice the density for peak calling and the 
      quaniles of the nonzero density values.
    }
  }
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
quants <- matrix(rnorm(12*5000), nrow = 5000, ncol = 12)
bounds <- cbind(mzmin = seq(from = 100, to = 1100, length.out = 5000),
                mzmax = seq(from = 100, to = 1100, length.out = 5000) + 0.1,
                scan.min = rep(10,5000), scan.max = rep(20, 5000),
                peaknum = 1:5000)
cmsobj <- new("CMSslice", assays = list(peakQuants = quants), 
                         rowData = DataFrame(bounds))
head(peakQuants(cmsobj))
}