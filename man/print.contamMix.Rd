\name{print.contamMix}
\alias{print.contamMix}
\title{
Print contamination estimates
}
\description{
This function displays a text summary of contamMix input data and MCMC chain
results (if runChain has been run).
}
\usage{
\method{print}{contamMix}(x, tabDelim=FALSE, ...)
}
\arguments{
  \item{x}{
    \code{contamMix} object (produced by \code{\link{loadSAM}} or
    \code{\link{runChains}})
  }
  \item{tabDelim}{
    If FALSE, output human-readable text.  If TRUE and \code{runChains}
    has been run, output a single line with the following tab-delimited
    columns:
    \enumerate{
      \item estimated error rate
      \item estimated maximum a posterior proportion authentic
      \item 2.5\% credible quantile for proportion authentic
      \item 97.5\% credible quantile for proportion authentic
      \item Gelman and Rubin diagnostic (\sQuote{potential scale
        reduction factor}) point estimate
      \item Gelman and Rubin diagnostic (\sQuote{potential scale
        reduction factor}) upper confidence limit
    }
  }
  \item{\dots}{
    (currently unused)
  }
}
\details{
Trace and posterior plots use MCMC samples only from the first chain,
regardless of how many chains were run.  Combining samples from multiple
chains is not technically legitimate.
}
\references{
As of May 2013, the best reference is the supplementary information of
Curr Biol. 2013 Apr 8;23(7):553-9.
}
\author{
  Philip Johnson <plfjohnson@emory.edu>
}

\seealso{
  \code{\link{contamMix}}
}
\examples{
\dontrun{
data = loadSAM(samFn="reads.sam", malnFn="aln.fa")
print(data)
res = runChains(data=data, nChains=3, nIter=10000, alpha=0.1)
print(res)
print(res, tabDelim=TRUE)
}
}
\keyword{ print }