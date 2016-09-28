\name{CMSproc-class}
\docType{class}
\alias{CMSproc-class}
\alias{densityEstimate}
\alias{densityCutoff}
\alias{densityQuantiles}
\alias{peakBounds}
\alias{peakQuants}
\alias{show,CMSproc-method}
\title{A class to hold chromatography-mass spectrometry preprocessing information.}
\description{
This class transforms information from raw mass spectral files to a matrix
of peak quantifications suitable for differential analysis. It first 
performs background correction and retention time alignment of the raw 
data. Next it performs weighted bivariate kernel density estimation to 
detect peaks.
}

\section{Slots}{
  \describe{
    \item{\code{colData}:}{
      a \code{DataFrame} of phenotype and sample information.
    }
    \item{\code{rawDT}:}{
      a \code{data.table} of raw spectral information.
    }
    \item{\code{mzParams}:}{
      a list containing the minimum and maximum M/Z value and number of 
      scans in each sample.
    }
    \item{\code{rtAlign}:}{
      a logical indicating whether the data has been retention time
      aligned or not.
    }
    \item{\code{bgcorrDT}:}{
      a \code{data.table} of background-corrected spectral information.
    }
    \item{\code{density}:}{
      a matrix with rows corresponding to M/Z values and columns 
      corresponding to scans containing the kernel density estimate.
    }
    \item{\code{densityCutoff}:}{
      the density cutoff value used to determine peak regions.
    }
    \item{\code{densityQuantiles}:}{
      a numeric vector containing the 100 percent quantiles of the
      nonzero density values.
    }
    \item{\code{peakBounds}:}{
      a matrix of information about the detected peaks : M/Z bounds,
      scan bounds, and ID number (\dQuote{peaknum}).
    }
    \item{\code{peakQuants}:}{
      a matrix of quantifications where rows correspond to peaks and 
      columns correspond to samples.
    }
  }
}

\section{Utility functions}{
  We have the following utility functions:

  \describe{
    \item{\code{show}:}{
      The show method; prints the object.
    }
    \item{\code{getEICS}:}{
      Gets extracted ion chromatograms (EICs) for the supplied M/Z ranges.
    }
    \item{\code{updatePeaks}:}{
      Given a new threshold value for the density estimate, repeats peak 
      detection and quantification and provides an updated \code{CMSproc} object.
    }
  }
}

\section{Accessors}{
  We have the following accessor functions:

  \describe{
    \item{\code{colData}:}{
      Gets the \code{DataFrame} containing phenotype and sample information.
    }
    \item{\code{densityEstimate}:}{
      Gets the matrix containing the density estimate.
    }
    \item{\code{densityCutoff}:}{
      Gets the value used to threshold the density for peak calling.
    }
    \item{\code{densityQuantiles}:}{
      Gets the quantiles of the nonzero values in the density estimate.
    }
    \item{\code{peakBounds}:}{
      Gets the M/Z bounds, scan bounds, and ID numbers for detected peaks.
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
cmsobj <- new("CMSproc", peakQuants = quants, peakBounds = bounds)
head(peakQuants(cmsobj))
}