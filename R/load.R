loadSAM <- function(samFn, malnFn, endogId = NULL, baseqThreshold = 30,
                    transverOnly = FALSE, trimBases = 0, saveMN = FALSE) {
  ## find our perl script
  sam2mn = system.file("exec/sam2mn.pl", package="contamMix")
  if (sam2mn == "") stop("could not find 'sam2mn.pl' script!")

  ## error checks
  if (is.null(samFn)  ||  !file.exists(samFn))
    stop("Cannot find SAM file '", samFn, "'.")
  if (is.null(malnFn)  ||  !file.exists(as.character(malnFn)))
    stop("Cannot find multiple alignment file '", malnFn, "'.")
  
  ## load data
  optnotnull <- function(opt,val) (
    if (is.null(val)  ||  val == ""  ||  val == FALSE)
    "" else paste(opt,val))
  cmd = paste("perl -w", sam2mn, "-bam", samFn, "-fa", malnFn,
      optnotnull("-ref", endogId), optnotnull("-baseq", baseqThreshold),
      optnotnull("-transverOnly", transverOnly),
      optnotnull("-trimBases", trimBases))
  if (!saveMN) {
      x = read.table(pipe(cmd), comment.char='', fill=TRUE)
  } else {
      mnFilename = paste(samFn,".mn", sep='');
      system(paste(cmd, ">", mnFilename))
      x = read.table(mnFilename, comment.char='', fill=TRUE)
  }

  data=x_loadMNtable(x)
  data$samFn=samFn
  data
}

x_loadMNtable <- function(x) {
  ## parse out error estimate
  if (x[nrow(x),1] != "#error:") { # independently estimated via other means
    stop("error parsing input data")
  }
  priorErrorRate = as.numeric(as.character(x[nrow(x),2]))
  if (priorErrorRate == 0) { priorErrorRate = 1e-6 } #avoid log(0)s
  if (priorErrorRate > 0.15) {
    stop("estimated error rate is implausibly high -- check for error in alignments")
  }

  ## reshape data into tensor (3-d matrix)
  numMT = (ncol(x)-1)/2
  mnMatrix = as.matrix(x[-nrow(x),-1])
  dim(mnMatrix) = c(nrow(mnMatrix), 2, numMT)
  dimnames(mnMatrix) = list(x[-nrow(x),1], c("M","N"), 0:(numMT-1))

  structure(list(numMT=numMT, mnMatrix=mnMatrix, e=priorErrorRate),
            class = "contamMix")
}
