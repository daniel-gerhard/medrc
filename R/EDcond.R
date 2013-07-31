EDcond <- function (object, respLev, cranef=NULL, interval = c("none", "delta", "fls", "tfls"), level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), type = c("relative", "absolute"), lref, uref, bound = TRUE, display = TRUE, logBase = NULL, ...){
  interval <- match.arg(interval)
  reference <- match.arg(reference)
  type <- match.arg(type)
  if ((type == "relative") && (bound)) {
    if (any(respLev <= 0 | respLev >= 100)) {
      stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
    }
  }
  EDlist <- object$fct$edfct
  if (is.null(EDlist)) {
    stop("ED values cannot be calculated")
  }
  indexMat <- object$indexMat
  parmMat <- object$parmMat
  
  pnames <- rownames(parmMat)
  
  #### cranef=NULL
  if (is.null(cranef)){
    randinfo <- lapply(object$fit$plist, function(x) x$random)
    ran <- ranef(object$fit)  
    
    # single vc
    rn <- names(ran)
    rnspl <- strsplit(rn, ".", fixed=TRUE)
    nrn <- sapply(rnspl, function(x) x[1])
    rest <- sapply(rnspl, function(x) x[2])
  } else {
    ran <- cranef
    if (is.null(rownames(ran))) rownames(ran) <- 1:nrow(ran)
    if (is.null(colnames(ran))) stop("cranef needs column names!")
    if (any(!(colnames(ran) %in% pnames))) stop("column names of cranef need to be a parameter name!")
    nrn <- colnames(ran)
  }
    
  mr <- matrix(0, ncol=length(pnames), nrow=nrow(ran), dimnames=list(rownames(ran), pnames))
  wn <- sapply(1:length(nrn), function(i) which(pnames == nrn[i]))
  mr[,wn] <- as.matrix(ran)
  
  # extend parmMat
  parmMat <- matrix(as.vector(apply(mr, 1, function(x) parmMat + x)), nrow=nrow(parmMat), dimnames=list(pnames, paste(rep(colnames(parmMat), times=nrow(mr)), rep(rownames(ran), each=ncol(indexMat)), sep=":")))
  
  indexMat <- matrix(1:length(parmMat), nrow=nrow(parmMat), dimnames=list(pnames, colnames(parmMat)))
  
  strParm0 <- sort(colnames(parmMat))
  curveNames <- colnames(parmMat)
  options(warn = -1)
  if (any(is.na(as.numeric(curveNames)))) {
    curveOrder <- order(curveNames)
  } else {
    curveOrder <- 1:length(curveNames)
  }
  options(warn = 0)
  strParm0 <- curveNames[curveOrder]
  indexMat <- indexMat[, curveOrder, drop = FALSE]
  parmMat <- parmMat[, curveOrder, drop = FALSE]
  strParm <- strParm0
  vcMat <- vcov(object)
  vcMat <- kronecker(diag(nrow(mr)), vcMat)
  ncolIM <- ncol(indexMat)
  indexVec <- 1:ncolIM
  lenPV <- length(respLev)
  noRows <- ncolIM * lenPV
  dimNames <- rep("", noRows)
  EDmat <- matrix(0, noRows, 2)
  oriMat <- matrix(0, noRows, 2)
  if (identical(length(unique(strParm)), 1)) {
    strParm[indexVec] <- rep("", ncolIM)
  } else {
    strParm <- paste(strParm, ":", sep = "")
  }
  rowIndex <- 1
  for (i in indexVec) {
    parmChosen <- parmMat[, i]
    parmInd <- indexMat[, i]
    varCov <- vcMat[parmInd, parmInd]
    for (j in 1:lenPV) {
      EDeval <- EDlist(parmChosen, respLev[j], reference = reference, type = type, ...)
      EDval <- EDeval[[1]]
      dEDval <- EDeval[[2]]
      oriMat[rowIndex, 1] <- EDval
      oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
      if (!is.null(logBase)) {
        EDval <- logBase^(EDval)
        dEDval <- EDval * log(logBase) * dEDval
      }
      EDmat[rowIndex, 1] <- EDval
      EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
      dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
      rowIndex <- rowIndex + 1
    }    
  }
  colNames <- c("Estimate", "Std. Error")
  if (interval == "delta") {
    intMat <- drc:::confint.basic(EDmat, level, object$type, df.residual(object), FALSE)
    intLabel <- "Delta method"
  }
  if (interval == "tfls") {
    intMat <- exp(drc:::confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2]/oriMat[, 1]), ncol = 2), level, object$type, df.residual(object), FALSE))
    intLabel <- "To and from log scale"
  }
  if (interval == "fls") {
    if (is.null(logBase)) {
      logBase <- exp(1)
      EDmat[, 1] <- exp(EDmat[, 1])
    }
    intMat <- logBase^(drc:::confint.basic(oriMat, level, object$type, df.residual(object), FALSE))
    intLabel <- "Back-transformed from log scale"
    EDmat <- EDmat[, -2, drop = FALSE]
    colNames <- colNames[-2]
  }
  if (identical(interval, "none")) {
    intLabel <- NULL
  } else {
    EDmat <- as.matrix(cbind(EDmat, intMat))
    colNames <- c(colNames, "Lower", "Upper")
  }
  dimnames(EDmat) <- list(dimNames, colNames)
  drc:::resPrint(EDmat, "Estimated effective doses", interval, intLabel, display = display)
}