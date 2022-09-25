\name{plot.contamMix}
\alias{plot.contamMix}
\title{
Plot diagnostics for contamination estimation
}
\description{
This function generates three types of plots used for assessing
MCMC convergence and quality of posterior density estimates.
}
\usage{
\method{plot}{contamMix}(x, which=1:3, ...)
}
\arguments{
  \item{x}{
    \code{contamMix} object (produced by \code{\link{runChains}})
  }
  \item{which}{
    Specifies which of the three plots to create (by default, all):
    \enumerate{
      \item Gelman diagonostic (see \code{\link{gelman.plot}})
      \item Trace of Pr{authentic} as a function of the MCMC iteration
      in the first chain only.
      \item Posterior kernel density estimate for Pr{authentic} for the
      first chain only.
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
res = runChains(data=data, nChains=3, nIter=10000, alpha=0.1)
plot(res, which=3)
}
}
\keyword{ hplot }
