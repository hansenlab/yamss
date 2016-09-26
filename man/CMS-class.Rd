\name{CMS-class}
\docType{class}
\alias{CMS-class}
\alias{getFileNames}
\alias{getRawPeakInfo}
\alias{densityEstimate}
\alias{densityCutoff}
\alias{densityQuantiles}
\alias{peakBounds}
\alias{peakQuants}
\alias{show,CMS-method}
\title{A class to hold mass spectrometry preprocessing information.}
\description{
This class transforms information from raw mass spectral files to a matrix
of peak quantifications suitable for differential analysis. It first 
performs background correction and retention time alignment of the raw 
data. Next it performs weighted bivariate kernel density estimation to 
detect peaks. Finally peak integration and differential analysis are
performed.
}

\section{Slots}{
  \describe{
    \item{\code{fileNames}:}{
      a character vector of filepaths to the raw data.
    }
    \item{\code{rtalign}:}{
      logical - was retention time alignment performed?
    }
    \item{\code{rawPeakInfo}:}{
      a list of raw spectral information for each scan in each sample.
    }
    \item{\code{rawPeakDT}:}{
      a \code{data.table} of raw spectral information.
    }
    \item{\code{bgcorrDT}:}{
      a \code{data.table} of background-corrected spectral information.
    }
    \item{\code{mzParams}:}{
      a list containing the minimum and maximum M/Z value and number of 
      scans in each sample.
    }
    \item{\code{bgSmooths}:}{
      a list containing M/Z region-specific smoothed background trends as 
      a function of scan for each sample.
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
    \item{\code{alignments}:}{
      a list of shifts for each sample in several M/Z regions that 
      maximizes the correlation between the extracted ion chromatograms 
      in those regions.
    }
    \item{\code{xicsRaw}:}{
      a list of extracted ion chromatograms corresponding to each peak 
      for each sample.
    }
    \item{\code{xicsImputed}:}{
      a list of interpolated XICs used for quantification. This is only 
      computed when retention time alignment is not performed.
    }
    \item{\code{peakBounds}:}{
      a matrix of information about the detected peaks : M/Z bounds,
      scan bounds, and ID number (\dQuote{blobnum}).
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
      detection and quantification and provides an updated \code{CMS} object.
    }
  }
}

\section{Accessors}{
  We have the following accessor functions:

  \describe{
    \item{\code{getFileNames}:}{
      gets filenames of raw data.
    }
    \item{\code{densityEstimate}:}{
      gets the matrix containing the density estimate.
    }
    \item{\code{densityCutoff}:}{
      gets the value used to threshold the density for peak calling.
    }
    \item{\code{densityQuantiles}:}{
      gets the quantiles of the nonzero values in the density estimate.
    }
    \item{\code{peakBounds}:}{
      gets the M/Z bounds, scan bounds, and ID numbers for detected peaks.
    }
    \item{\code{peakQuants}:}{
      gets the matrix of peak quantifications (rows: peaks, columns: samples).
    }
  }
}

\examples{
quants <- matrix(rnorm(12*5000), nrow = 5000, ncol = 12)
bounds <- cbind(mzmin = seq(from = 100, to = 1100, length.out = 5000),
                mzmax = seq(from = 100, to = 1100, length.out = 5000) + 0.1,
                scan.min = rep(10,5000), scan.max = rep(20, 5000),
                blobnum = 1:5000)
cmsobj <- new("CMS", peakQuants = quants, peakBounds = bounds)
head(peakQuants(cmsobj))
}