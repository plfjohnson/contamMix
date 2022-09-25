## ---------------------------------------------------------------------------
## pre : x and y coords intended for plotting, desired output dpi
## post: which coords should be plotted for optimal display at specified dpi
x_thinXY <- function (x, y=NULL, type='p', dpi=20, xlim=NULL, ylim=NULL) {
  xy = xy.coords(x, y)
  if (is.null(xlim)) xlim = range(xy$x)
  if (is.null(ylim)) ylim = range(xy$y)
  uin = c(diff(xlim), diff(ylim))/par("pin") #usr coord / inch
  udt = uin / dpi #usr coord / dot
  dtx = round(xy$x / udt[1])
  dty = round(xy$y / udt[2])
  if (type == 'p') {## if plotting points, just eliminate all dups
    w=!duplicated(cbind(dtx,dty))
  } else if (type == 'l') { ## if plotting lines, eliminate all *pairs* of dups
    n = length(dtx)
    z=cbind(dtx[-n],dty[-n],dtx[-1],dty[-1])
    ## line is same either direction
    rv = z[,1] > z[,3]  |  (z[,1] == z[,3]  &  z[,2] > z[,4])
    z[rv,] = z[rv,c(3,4,1,2)]
    gdPairs = !duplicated(z)  &  (dtx[-n] != dtx[-1]  |  dty[-n] != dty[-1])
    w=c(gdPairs, FALSE)  |  c(FALSE, gdPairs)
  }
  xy$x = xy$x[w]
  xy$y = xy$y[w]
  return(xy)
}

x_AnalyzeChains <- function(chains, focalChain=1) {
  ## crude identification of burn-in
  ## could just use "lik", but empirical results suggest "0" is often better..
  burnIn = max(sapply(chains, function(x) (max(
    which(x[,"0"]   > median(x[(nrow(x) %/% 2):nrow(x),"0"])  )[1],
    which(x[,"lik"] > median(x[(nrow(x) %/% 2):nrow(x),"lik"]))[1]))
    ))

  if (is.na(burnIn)) {
    warning("MCMC not run long enough to guestimate burn in; proceeding using no burn in for diagonostic purposes.")
    burnIn=1;
  }

  ## arbitrarily designate one run to analyze
  chain = chains[[focalChain]]
  samples = chain[-(1:burnIn),"0"]

  ## estimate posterior density after burn-in
  d=density(samples, cut=0);
  mcmcObj = coda::mcmc.list(lapply(chains, function(x) (coda::mcmc(x[,"0"]))))

  list(burnIn=burnIn, chain=chain, samples=samples, dens=d,
       MAP = d$x[which.max(d$y)], ci=quantile(samples, c(0.025, 0.975)),
       mcmcObj=mcmcObj)
}

## ---------------------------------------------------------------------------
## Public functions below

print.contamMix <- function(x, tabDelim = FALSE, ...) {
  if (!tabDelim) {
    cat("Data from:\n  ", x$samFn, "\nconsist of", nrow(x$mnMatrix), "reads,",
        "considered as a mixture of", x$numMT, "possible genomes.\n")

    contamEvidence = apply(x$mnMatrix, 1, function(d) (any(d["M",1]<d["M",-1])))
    cat("Crude contamination upper bound: ", sum(contamEvidence), " out of ",
        nrow(x$mnMatrix), " reads (",
        signif(100*sum(contamEvidence)/nrow(x$mnMatrix), 2),
        "%) match another\n",
        "genome better than the consensus (error rate is ",
        x$e, ").\n\n", sep='')

    if (is.null(x$chains)) {
      cat("[mcmc estimate not available since chain hasn't been run]\n\n")
    }
  }
  
  if (is.null(x$chains)) return(invisible())

  ## arbitrarly focus on one chain
  anal = x_AnalyzeChains(x$chains, focalChain=1)

  if (tabDelim) {
    cat(x$e, anal$MAP, anal$ci,
        if (length(x$chains) > 1) coda::gelman.diag(anal$mcmcObj)$psrf else "",
        sep="\t");
    cat("\n")
  } else {
    cat("quantiles from n = ", length(anal$samples),
        " samples (after discarding burnin of ", anal$burnIn, ")\n", sep='')
    print(anal$ci)
    cat("MAP authentic:", anal$MAP, "\n");
    if (length(x$chains) > 1) print(coda::gelman.diag(anal$mcmcObj))
  }
  
  invisible()
}


plot.contamMix <- function(x, which=1:3, ...) {
  if (is.null(x$chains)) stop("MCMC not yet run (use runChain)")
  
  anal = x_AnalyzeChains(x$chains, focalChain=1)
  
  ##convergence diagnostic, if applicable
  if (1 %in% which  &&  length(x$chains) > 1) {
    coda::gelman.plot(anal$mcmcObj, auto.layout=FALSE, ylim=c(1,2))
  }

  if (2 %in% which) {
    plot(x_thinXY(anal$chain[,"0"], dpi=60), type='l', xlab='Iteration', ylab='Pr(authentic)')
    abline(v=anal$burnIn, col = "blue", lty=3)
  }

  if (3 %in% which) {
    plot(anal$dens, type='n', main=NA, xlab="Pr(authentic)", xlim=c(0,1))
    ciRegion = anal$dens$x >= anal$ci[1] & anal$dens$x <= anal$ci[2]
    polygon(cbind(c(anal$ci[1], anal$dens$x[ciRegion], anal$ci[2]),
                  c(0, anal$dens$y[ciRegion], 0)), border=NA, col='gray')
    lines(anal$dens) #plot density on top of shaded region
    abline(v=anal$MAP, lty=2, col='red')
    text(x=anal$MAP, y=mean(par()$usr[3:4]), labels=
         paste("~", format(anal$MAP, nsmall=2, digits=2), sep=''),
         col='red', cex=1.2, adj=c(1.2, .5), xpd=NA)
  }
  
  invisible()
}
