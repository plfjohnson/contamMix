x_gibbs.step <- function(D, prop, alpha) {
  .Call("gibbs_step", D, prop, alpha)
}

## rdirichlet code posted by Ben Bolker to R-News on Fri Dec 15
## 2000. See http://tolstoy.newcastle.edu.au/R/help/00b/2175.html. Ben
## attributed the code to Ian Wilson i.wilson@maths.abdn.ac.uk.
x_rdirichlet <- function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

x_runChain <- function(y, nIter, alpha, recordAll) {
  res = matrix(0, nrow=nIter, ncol=1 + ifelse(recordAll, ncol(y),1))
  colnames(res) = c("lik", 0:(ncol(res)-2));
  ## initialize chain at random proportions
  mcmcProp = as.vector(x_rdirichlet(1, rep(alpha, ncol(y))))
  ## iterate MCMC
  for (i in 1:nIter) {
    res[i,]=c(sum(log(y %*% mcmcProp)),
         if (recordAll) mcmcProp else mcmcProp[1])
    mcmcProp = x_gibbs.step(y, mcmcProp, alpha)
  }
  return(res)
}

runChains <- function(data, nChains, nIter, alpha = 0.1, recordAll = FALSE) {
  if (class(data) != "contamMix") stop("data is of wrong type")
  if (nChains < 1) stop("nChains must be a positive number")
  if (nIter < 1) stop("nIter must be a positive number")
  if (alpha < 0) stop("alpha must be non-negative number")

  ## incorporate errors into data
  y = (1-data$e)^data$mnMatrix[,"M",] * data$e^data$mnMatrix[,"N",]
  if (is.vector(y)) { y = t(y) } #edge case if only one read

  ## run multiple chains, ideally in parallel
  data$chains =
    parallel::mclapply(1:nChains, function(i) (
      x_runChain(y, nIter, alpha, recordAll) ))
  data
}
