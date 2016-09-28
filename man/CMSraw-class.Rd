\name{CMSraw-class}
\docType{class}
\alias{CMSraw-class}
\alias{show,CMSraw-method}
\alias{colData,CMSraw-method}
\title{A class to hold chromatography-mass spectrometry raw data and metadata.}
\description{
This class saves the data from raw mass spectral data files in a 
\code{data.table} and is used in subsequent processing. Objects of 
this class are created by the \code{readMSdata} function.
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
  }
}

\section{Accessors}{
  We have the following accessor functions:

  \describe{
    \item{\code{colData}:}{
      Gets the \code{DataFrame} containing phenotype and sample information.
    }
  }
}

\examples{
if (require(mtbls2)) {
data(mtbls2)
filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1", recursive = TRUE, full.names = TRUE)[1:2]
colData <- DataFrame(group = rep(c("wild-type", "mutant"), each = 4))[1:2,]
cmsRaw <- readMSdata(files = files, colData = colData, verbose = TRUE)
colData(cmsRaw)
}
}