\name{loadSAM}
\alias{loadSAM}
\title{
Process SAM or BAM file.
}
\description{
Loads MT data from SAM or BAM file and combines with information in
multiple alignment of endogenous sequence + potential contaminants.
Requires PERL to be installed and in your path.
}
\usage{
loadSAM(samFn, malnFn, endogId = NULL, baseqThreshold = 30, transverOnly =
        FALSE, trimBases = 0, saveMN = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samFn}{
    path to SAM or BAM file.  If BAM, also requires \code{samtools} to be
    installed and in path.  Reads \emph{must} be aligned to endogenous
    sequence (presumably the consensus).
  }
  \item{malnFn}{
    path to FASTA-formatted multiple alignment file containing
    endogenous genome + potential contaminants.
  }
  \item{endogId}{
    If NULL, assumes the endogenous ID in the alignment file exactly
    matches the reference ID in the SAM file.  If these do not match
    (but are, in fact, the same sequence), then use \code{endogId} to
    specify the name in the alignment file.
  }
  \item{baseqThreshold}{
    Discards any base with quality score lower than this value.
  }
  \item{transverOnly}{
    If TRUE, only use sites in which the alleles form a transversion --
    this loses power due to the reduction in number of sites, but avoids
    any potential biases due to aDNA damage.
  }
  \item{trimBases}{
    Trim this many bases from both ends of each read, in an attempt to
    reduce the effect of aDNA damage without reducing the number of
    informative sites (i.e., the \code{transverOnly} option).
  }
  \item{saveMN}{
    Whether to save MN intermediate file (otherwise just uses pipe so never
    explicitly saves to disk).  If true, will save to filename
    \sQuote{<samFn>.mn}.
  }
}
\details{
  Under the hood, this function runs the perl script \code{sam2mn.pl}, which is
  included with this R package.  This script collapses the input
  alignment data into simple summary statistics: for each read, the #
  matches and # nonmatches to each genome.  For a dataset with many
  reads, this script may take a while.

  The returned object can be \code{print}ed for a human-readable summary
  of the data before handing off to \code{runChains}.  Or the object could
  be \code{save}d to a file for future use.
}
\value{
  \code{contamMix} object, which consists of
  \item{numMT}{The number of mitochondrial genomes in the multiple
    alignment}
  \item{mnMatrix}{3-d matrix of match/nonmatch counts for each read and
    genome}
  \item{e}{Estimated error rate}
}
\references{
As of May 2013, the best reference is the supplementary information of
Curr Biol. 2013 Apr 8;23(7):553-9.
}
\author{
  Philip Johnson <plfjohnson@emory.edu>
}

\seealso{
  \code{\link{contamMix}} and specifically \code{\link{runChains}}
}
\examples{
\dontrun{
data = loadSAM(samFn="reads.sam", malnFn="aln.fa")
print(data)
}
}
\keyword{ IO }

