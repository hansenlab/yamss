\name{cms-class}
\docType{class}
\alias{cms-class}
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
    \item{\code{files}:}{
       a character vector of filepaths to the raw data
    }
    \item{\code{classes}:}{
       a vector of sample class information corresponding to \code{files}
    }
    \item{\code{rawpeakinfo}:}{
       a list of raw spectral information for each scan in each sample
    }
    \item{\code{mzParams}:}{
       a list containing the minimum and maximum M/Z value and number of 
       scans in each sample
    }
    \item{\code{bgsmooths}:}{
       a list containing M/Z region-specific smoothed background trends as 
       a function of scan for each sample
    }
    \item{\code{dens}:}{
       a matrix with rows corresponding to M/Z values and columns 
       corresponding to scans containing the kernel density estimate
    }
    \item{\code{dcutoff}:}{
       the density cutoff value used to determine peak regions
    }
    \item{\code{densquants}:}{
       a numeric vector containing the quantiles of the nonzero density 
       values
    }
    \item{\code{alignments}:}{
       a list of shifts for each sample in several M/Z regions that 
       maximizes the correlation between the extracted ion chromatograms 
       in those regions
    }
    \item{\code{blobs}:}{
       a matrix of information about the peaks (high density "blobs") 
       detected: M/Z bounds, scan bounds, and ID number ("blobnum")
    }
    \item{\code{xicsRaw}:}{
       a list of extracted ion chromatograms corresponding to each peak 
       for each sample
    }
    \item{\code{xicsImputed}:}{
       a list of interpolated XICs used for quantification. This is only 
       computed when retention time alignment isn't performed
    }
    \item{\code{quants}:}{
       a matrix of quantifications where rows correspond to peaks and 
       columns correspond to samples
    }
    \item{\code{diffrep}:}{
       a \code{data.frame} containing the results of a differential 
       analysis
    }
  }
}