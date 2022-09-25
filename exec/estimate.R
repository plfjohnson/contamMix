library(contamMix)
library(getopt)
options(mc.cores = parallel::detectCores())


## load arguments from command line (if any) & set defaults
opt = getopt(rbind(
  c('samFn', NA, 1, "character"),
  c('malnFn', NA, 1, "character"),
  c('consId', NA, 1, "character"),
  c('figure', NA, 1, "character"), #if !interactive, save figure to this file
  c('nIter', NA, 1, "integer"),
  c('nChains', NA, 1, "integer"),
  c('alpha', NA, 1, "double"),
  c('baseq', NA, 1, "integer"),
  c('transverOnly', NA, 0, "logical"),
  c('trimBases', NA, 1, "integer"),
  c('saveData', NA, 1, "character"),#save MCMC data for manual diag
  c('saveMN', NA, 0, "logical"),    #save MN intermediate for debugging
  c('recordAll', NA, 0, "logical"), #record *all* proportions from MCMC chain
  c('tabOutput', NA, 0, "logical")) #simple tab-delimited output
  )
if (is.null(opt$nIter))  opt$nIter = 50000
if (is.null(opt$nChains)) opt$nChains = 3
if (is.null(opt$alpha)) opt$alpha = 0.1
if (is.null(opt$baseq)) opt$baseq = 30; if (opt$baseq == 0) opt$baseq = NULL
if (is.null(opt$recordAll)) opt$recordAll = FALSE
if (is.null(opt$saveMN)) opt$saveMN = FALSE
if (is.null(opt$tabOutput)) opt$tabOutput = FALSE

if (!interactive()  &&  (is.null(opt$samFn)  ||  is.null(opt$malnFn))) {
  cat("Usage: ./estimate.R --samFn <data.sam>  --malnFn <alignment.fa> [--nIter <50000>] [--nChains <3>] [--alpha <0.1>] [--figure <outputfigure.pdf>] [--baseq <30>] [--transverOnly] [--trimBases <0>] [--saveMN] [--tabOutput] ...\n\n",
      "Required parameter:\n",
      "\t--samFn    --> SAM/BAM data file aligned to consensus\n",
      "\t--malnFn   --> FASTA-format multiple alignment of consensus and potential contaminants\n",
      "Optional parameters:\n",
      "\t--consId --> consensus id in multiple alignment (if not supplied, assumes this is identical to the reference id in the SAM/BAM)\n",
      "\t--nIter --> number of iterations in Markov chain (default 50000)\n",
      "\t--nChains --> number of MCMC chains to run from different random starting parameters, which are used for Gelman-Rubin convergence testing (default 3; uses multiple processors if available)\n",
      "\t--alpha --> hyperparameter to Dirichlet prior distribution; may need to be tweaked if MCMC is mixing poorly (default 0.1)\n",
      "\t--figure --> if supplied, generates a PDF figure with 3 panels: convergence of gelman diagnostic (if --nChains>1),  Pr(authentic) as a function of MC iteration, estimated posterior density for Pr(authentic)\n",
      "\t--baseq --> base quality threshold below which to discard data (default 30; 0 signals NO threshold)\n",
      "\t--transverOnly --> only use sites with transversions (avoids potential for bias from aDNA damage, but has significantly less power)\n",
      "\t--trimBases --> trim this # of bases from ends of sequence (reducing effect of aDNA damage)\n",
      "\t--saveData --> save chain data to specified file (in .Rdata format) for manual diagnostics\n",
      "\t--saveMN --> save MN intermediate file for manual debugging (will use filename '<samFn>.mn')\n",
      "\t--tabOutput --> output a single line of text with the following tab-separated values: <inferred-error-rate> <MAP-authentic> <2.5% authentic> <97.5% authentic> <gelman diagnostic> <gelman diag upper bound>\n\n")
  cat("NOTE: Not converged if Gelman diagnostic >1.1.  Initial runs on any new data should request a --figure to visually check that the Markov chain is well-behaved.\n\n")
  q()
}

## test only
#data = loadSAM(samFn="~/data/emh/mt_only/Paglicci23_sequence_merged_quality_Monly-human_MT.sort.q30_l35_uniq.bam", malnFn="~/data/emh/mt_only/aln.mafft", opt$alnFn, endogId="Paglicci23"

## do the work
data = loadSAM(samFn=opt$samFn, malnFn=opt$malnFn, endogId=opt$consId, baseqThreshold=opt$baseq, transverOnly=opt$transverOnly, trimBases=opt$trimBases, saveMN=opt$saveMN)
res = runChains(nChains=opt$nChains, nIter=opt$nIter, data=data, alpha=opt$alpha)

## results in text format
print(res, tabDelim=opt$tabOutput)

## binary save of results
if (!is.null(opt$saveData)) {
  if (!grepl(".Rdata$", opt$saveData))
      opt$saveData = paste(opt$saveData, ".Rdata", sep="")
  if (!opt$tabOutput) cat("Saving data & results to", opt$saveData, "\n")
  save(res, file=opt$saveData)
}

## graphical results
if (!is.null(opt$figure)) {
  if (!grepl(".pdf$", opt$figure)) opt$figure = paste(opt$figure, ".pdf", sep="")
  if (!opt$tabOutput) cat("Saving figure to", opt$figure, "\n")
  
  dev.new(file=opt$figure, width=4, height=6, pointsize=10)
  par(mfrow=c(3,1), cex=1, mex=0.5, mar=c(5,5,1,1), oma=c(0,0,3,0))
  plot(res, which=1:3)
  mtext(opt$samFn, side=3, line=1, outer=TRUE, font=2)
  mtext(bquote(paste("(", epsilon, " = ", .(signif(res$e,2)), ")")), side=3, line=-1, outer=TRUE)
  invisible(dev.off())
}
