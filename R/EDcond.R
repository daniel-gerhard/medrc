EDcond <- function (object, respLev, rlevel=Q, interval = c("none", "delta", "fls", "tfls"), level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), type = c("relative", "absolute"), lref, uref, bound = TRUE, display = TRUE, logBase = NULL, ...){
  Q <- object$fit$dims$Q
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
  mcoefs <- cbind(coef(object$fit, level=rlevel))
  indexMat <- matrix(1:prod(dim(mcoefs)), nrow=ncol(mcoefs), ncol=nrow(mcoefs))
  parmMat <- t(mcoefs)
  
  
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
  vcMat <- kronecker(diag(ncol(parmMat)), vcov(object, od = od, pool = pool))
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
  lenIV <- length(indexVec)
  dEDmat <- matrix(0, lenPV * lenIV, nrow(vcMat))
  for (i in indexVec) {
    parmChosen <- parmMat[, i]
    parmInd <- indexMat[, i]
    varCov <- vcMat[parmInd, parmInd]
    for (j in 1:lenPV) {
      EDeval <- EDlist(parmChosen, respLev[j], reference = reference, 
                       type = type)
      EDval <- EDeval[[1]]
      dEDval <- EDeval[[2]]
      dEDmat[(i - 1) * lenPV + j, parmInd] <- dEDval
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
    intMat <- eval(parse(text="drc:::confint.basic(EDmat, level, object$type, df.residual(object), FALSE)"))
    intLabel <- "Delta method"
  }
  if (interval == "tfls") {
    intMat <- eval(parse(text="exp(drc:::confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2]/oriMat[, 1]), ncol = 2), level, object$type, df.residual(object), FALSE))"))
    intLabel <- "To and from log scale"
  }
  if (interval == "fls") {
    if (is.null(logBase)) {
      logBase <- exp(1)
      EDmat[, 1] <- exp(EDmat[, 1])
    }
    intMat <- eval(parse(text="logBase^(drc:::confint.basic(oriMat, level, object$type, df.residual(object), FALSE))"))
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
  eval(parse(text="drc:::resPrint(EDmat, 'Estimated effective doses', interval, intLabel, display = display)"))
}
 