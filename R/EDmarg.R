EDmarg <- function (object, respLev, interval = c("none", "delta", "fls", "tfls"), level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), type = c("relative", "absolute"), nGQ=5, rfinterval=c(-1000, 1000), lref, uref, bound = TRUE, display = TRUE, logBase = NULL, ...){
  require(statmod)
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
  
  # construct integration grid
  pnames <- rownames(parmMat)
  vc <- VarCorr(object$fit)
  nv <- nrow(vc)
  if (nv == 2){
    nvc <- names(ranef(object$fit))
  } else{
    nvc <- rownames(vc)[-nv] 
  }
  ###
  # single vc per parameter
  rnspl <- strsplit(nvc, ".", fixed=TRUE)
  nrn <- sapply(rnspl, function(x) x[1])
  rest <- sapply(rnspl, function(x) x[2])
  
  std <- as.numeric(vc[-nv,2])
  cormat <- diag(length(std))
  if (ncol(vc) > 2){    
    cormat[upper.tri(cormat)] <- cormat[lower.tri(cormat)] <- na.omit(as.numeric(vc[-c(1,nv),-c(1,2)]))
  }
  gq <- gauss.quad.prob(n=nGQ, dist="normal", sigma=1)
  
  eg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$nodes", nv-1), collapse=","), ")", sep="")))
  weg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$weights", nv-1), collapse=","), ")", sep="")))
  w <- apply(weg, 1, function(x) prod(x))
  
  ee <- eigen(cormat)
  A <- ee$vectors %*% diag(sqrt(ee$values))
  z <- data.frame(t(std*(A %*% t(as.matrix(eg)))))
  
  intgrid <- matrix(0, ncol=length(pnames), nrow=nrow(z))
  wn <- sapply(1:length(nrn), function(i) which(pnames == nrn[i]))
  intgrid[,wn] <- as.matrix(z)
  colnames(intgrid) <- pnames
  ###   
  
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
  
  # ED estimation by root finding based on marginal prediction
  EDlistm <- function(parmChosen, respLev, reference=reference, type=type, intgrid=intgrid, intweights=intweights, rfinterval=rfinterval){
    parm <- object$fct$fixed
    parm[is.na(parm)] <- parmChosen
    p <- 100-drc:::EDhelper(parmChosen, respLev, reference = reference, type = type)
    tval <- parm[2] + (parm[3]-parm[2])*(p/100)
    mint <- function(d) sum(na.omit(w*apply(intgrid, 1, function(x) object$fct$fct(d, rbind(parmChosen + x)))))-tval 
    myenv <- new.env()
    assign("tval", tval, envir = myenv)
    assign("parmChosen", parmChosen, envir = myenv) 
    assign("rfinterval", rfinterval, envir = myenv) 
    ede <- suppressWarnings(numericDeriv(quote(uniroot(mint, interval=rfinterval)$root), "parmChosen", myenv))
    out <- list()
    out[[1]] <- ede
    out[[2]] <- as.vector(attr(ede, "gradient"))
    return(out)
  }
  
  rowIndex <- 1
  for (i in indexVec) {
    parmChosen <- parmMat[, i]
    parmInd <- indexMat[, i]
    varCov <- vcMat[parmInd, parmInd]
    for (j in 1:lenPV) {
      EDeval <- EDlistm(parmChosen, respLev[j], reference = reference, type = type, intgrid=intgrid, intweights=intweights, rfinterval=rfinterval)
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