\name{contamMix-package}
\alias{contamMix-package}
\alias{contamMix}
\docType{package}
\title{
Estimate contamination proportion in MT reads
}
\description{
Assumes contamination <50\% and models data as a mixture of endogenous
and potential contaminant mitochondrial genomes.  Estimates mixture
proportions uses MCMC with a Dirichlet prior.
}
\details{
\tabular{ll}{
Package: \tab contamMix\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-15\cr
License: \tab GPL-3\cr
}

Mitochondrial genome reads should be aligned to their consensus sequence
(which is approximately the true sequence, assuming <50\% contamination).  Then the consensus genome and potential contaminant genomes must be multiply-aligned to each other (e.g., using a program like muscle or mafft).

Given these two pieces of input data, first call \code{\link{loadSAM}}
to load the data.  Then call \code{\link{runChains}} to perform the MCMC
estimation.  Finally output the results using \code{\link{print.contamMix}} or
\code{\link{plot.contamMix}}.

Perhaps the easiest way to perform estimation is to run the
\sQuote{estimate.R} command line script included in the exec directory
of this package.
}
\author{
Philip Johnson <plfjohnson@emory.edu>
}
\references{
As of May 2013, the best reference is the supplementary information of
Curr Biol. 2013 Apr 8;23(7):553-9.
}
\seealso{
  \code{\link{loadSAM}}, \code{\link{runChains}},
  \code{\link{print.contamMix}}, \code{\link{plot.contamMix}}
}
\keyword{ package }
\examples{
\dontrun{
data = loadSAM(samFn="reads.sam", malnFn="aln.fa")
res = runChains(data=data, nChains=3, nIter=10000, alpha=0.1)
print(res)
plot(res, which=3)
}
} 
