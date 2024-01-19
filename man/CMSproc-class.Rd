\name{CMSproc-class}
\docType{class}
\alias{CMSproc-class}
\alias{densityEstimate}
\alias{densityQuantiles}
\alias{show,CMSproc-method}
\title{A class to hold chromatography-mass spectrometry preprocessing information.}
\description{
This class builds on the \code{CMSraw} class to additionally store 
background-corrected intensities as well as the bivariate kernel density 
estimate.
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
    \item{\code{densityQuantiles}:}{
      a numeric vector containing the 100 percent quantiles of the
      nonzero density values.
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
    \item{\code{plotDensityRegion}:}{
      Makes an image plot of the density estimate in a specified M/Z and 
      scan region.
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
    \item{\code{densityQuantiles}:}{
      Gets the quantiles of the nonzero values in the density estimate.
    }
  }
}

\examples{

## Construct a completely fake example
densmat <- matrix(rnorm(600), nrow = 20, ncol = 30)
colnames(densmat) <- 1:ncol(densmat)
rownames(densmat) <- seq(350, by = 0.005, length.out = nrow(densmat))
cmsobj <- new("CMSproc", density = densmat)
head(densityEstimate(cmsobj))

## Takes about 20s to run
\donttest{

data(cmsRawExample)
cmsProc <- bakedpi(cmsRawExample,
                    dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                    dortalign = TRUE, mzsubset = c(500,510))
cmsProc
}
}