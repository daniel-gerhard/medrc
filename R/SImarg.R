SImarg <- function(object, percVec, percMat = NULL, compMatch = NULL, reverse = FALSE, interval = c("none", "delta", "fieller", "fls"), level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), type = c("relative", "absolute"), nGQ=5, rfinterval=c(-1000, 1000), display = TRUE, logBase = NULL, ...){  
  interval <- match.arg(interval)
  reference <- match.arg(reference)
  type <- match.arg(type)
  if ((is.null(logBase)) && (interval == "fls")) {
    stop("Argument 'logBase' not specified for interval = 'fls'")
  }
  if ((type == "relative") && any(percVec <= 0 | percVec >= 100)) {
    stop("Percentages outside the interval [0, 100] not allowed")
  }
  if (missing(compMatch)){
    matchNames <- FALSE
  } else {
    matchNames <- TRUE
  }
  lenPV <- length(percVec)
  curveNames <- colnames(object$parmMat)
  
  # ED estimation by root finding based on marginal prediction
  edfct <- function(parmChosen, respLev, reference, type, intgrid, intweights, rfinterval){
    parm <- object$fct$fixed
    parm[is.na(parm)] <- parmChosen
    parm2 <- object$fct$fixed
    
    p <- 100-eval(parse(text="drc:::EDhelper(parm, respLev, reference = reference, type = type)"))
    cip <- if ("c" %in% colnames(intgrid)) intgrid[,"c"] else rep(0, nrow(intgrid))
    dip <- if ("d" %in% colnames(intgrid)) intgrid[,"d"] else rep(0, nrow(intgrid))    
    mint <- function(d) sum(na.omit(w * sapply(1:nrow(intgrid), function(x){      
      pc <- if (is.na(parm2[2])) parmChosen["c"] else parm[2]
      pd <- if (is.na(parm2[3])) parmChosen["d"] else parm[3]
      object$fct$fct(d, rbind(parmChosen + intgrid[x,])) - ((pc + cip[x]) + ((pd + dip[x]) - (pc + cip[x])) * (p/100))      
    }) ))
    myenv <- new.env()
    assign("parmChosen", parmChosen, envir = myenv) 
    assign("rfinterval", rfinterval, envir = myenv) 
    ede <- suppressWarnings(numericDeriv(quote(uniroot(mint, interval=rfinterval)$root), "parmChosen", myenv))
    out <- list()
    out[[1]] <- ede
    out[[2]] <- as.vector(attr(ede, "gradient"))
    return(out)
  }
  
  indexMat <- object$indexMat
  parmMat <- object$parmMat
  
  pnames <- rownames(parmMat)
  
  # Pinheiro function
  varRan <- function(object, level = 1){
    sigE <- object$sig^2
    sigE*pdMatrix(object$modelStruct$reStruct)[[level]]
  }
  
  # extract variance components
  ranefs <- ranef(object$fit)
  nvc <- length(object$fit$modelStruct$reStruct)
  vclist <- lapply(1:nvc, function(i) varRan(object$fit, level=i))
  strspl <- lapply(vclist, function(x) strsplit(names(diag(x)), ".", fixed=TRUE))
  nrnl <- lapply(strspl, function(x) sapply(x, function(x) x[1]))
  
  ###
  gq <- gauss.quad.prob(n=nGQ, dist="normal", sigma=1)  
  
  igwlist <- lapply(1:nvc, function(v){
    nrn <- nrnl[[v]]
    nv <- length(diag(vclist[[v]]))
    
    eg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$nodes", nv), collapse=","), ")", sep="")))
    weg <- eval(parse(text=paste("expand.grid(", paste(rep("gq$weights", nv), collapse=","), ")", sep="")))
    w <- apply(weg, 1, function(x) prod(x))
    
    cfvv <- chol(vclist[[v]])
    z <- as.matrix(eg) %*% cfvv
    
    intgrid <- matrix(0, ncol=length(pnames), nrow=nrow(z))
    wn <- sapply(1:length(nrn), function(i) which(pnames == nrn[i]))
    intgrid[,wn] <- as.matrix(z)
    colnames(intgrid) <- pnames
    return(list(intgrid, w))
  })
  ###   
  iglist <- lapply(igwlist, function(x) x[[1]])
  ni <- sapply(iglist, nrow)
  vllist <- lapply(1:length(ni), function(i) apply(iglist[[i]], 2, function(ir){
    if (i == 1) return(rep(ir, each=prod(ni[(i+1):length(ni)])))
    if (i > 1 & i < length(ni)) return(rep(rep(ir, each=prod(ni[(i+1):length(ni)])), times=prod(ni[1:(i-1)])))
    if (i == length(ni)) return(rep(ir, times=prod(ni[1:(i-1)])))    
  }))
  arr <- simplify2array(vllist)
  intgrid <- apply(arr, c(1,2), sum)
  
  wlist <- lapply(igwlist, function(x) x[[2]])
  wmat <- sapply(1:length(ni), function(i){
    if (i == 1) return(rep(wlist[[i]], each=prod(ni[(i+1):length(ni)])))
    if (i > 1 & i < length(ni)) return(rep(rep(wlist[[i]], each=prod(ni[(i+1):length(ni)])), times=prod(ni[1:(i-1)])))
    if (i == length(ni)) return(rep(wlist[[i]], times=prod(ni[1:(i-1)])))    
  })
  w <- apply(wmat, 1, prod)
  
  sifct <- createsifctm(edfct, logBase, identical(interval,"fls"), indexMat, length(coef(object)), intgrid, w, rfinterval)
  
  options(warn = -1)
  if (any(is.na(as.numeric(curveNames)))){
    curveOrder <- order(curveNames)
  } else {
    curveOrder <- 1:length(curveNames)
  }
  options(warn = 0)
  strParm0 <- curveNames[curveOrder]
  indexMat <- indexMat[, curveOrder, drop = FALSE]
  lenEB <- ncol(indexMat)
  #sifct <- createsifct(object$fct$edfct, logBase, identical(interval, "fls"), indexMat, length(coef(object)))
  parmMat <- parmMat[, curveOrder, drop = FALSE]
  strParm <- strParm0
  varMat <- vcov(object)
  numComp <- (lenPV * (lenPV - 1)/2) * (lenEB * (lenEB - 1)/2)
  matchVec <- rep(TRUE, numComp)
  rNames <- rep("", numComp)
  oriMat <- matrix(0, numComp, 2)
  degfree <- df.residual(object)
  rowIndex <- 1
  pairsMat <- combinations(lenEB, 2)
  if (is.null(percMat)) {
    percMat <- combinations(lenPV, 2)
  }
  if (reverse){
    pairsMat <- pairsMat[, 2:1, drop = FALSE]
    percMat <- percMat[, 2:1, drop = FALSE]
  }
  appFct1 <- function(percVal){
    apply(pairsMat, 1, siInnerm, pVec = percVec[percVal], compMatch = compMatch, object = object, indexMat = indexMat, parmMat = parmMat, varMat = varMat, level = level, reference = reference, type = type, sifct = sifct, interval = interval, degfree = degfree, logBase = logBase, intgrid=intgrid, intweights=w, rfinterval=rfinterval)
  }
  SImat <- matrix(apply(percMat, 1, appFct1), nrow = nrow(pairsMat) * nrow(percMat), byrow = TRUE)
  
  appFct2 <- function(percVal) {
    apply(pairsMat, 1, function(indPair, percVal){
      paste(strParm[indPair[1]], "/", strParm[indPair[2]], ":", percVec[percVal[1]], "/", percVec[percVal[2]], sep = "")
    }, percVal = percVal)
  }
  rownames(SImat) <- apply(percMat, 1, appFct2)
  appFct3 <- function(percVal){
    apply(pairsMat, 1, function(indPair, percVal){
      (is.null(compMatch) || all(c(strParm[indPair[1]], strParm[indPair[2]]) %in% compMatch))
    })
  }
  SImat <- SImat[as.vector(apply(percMat, 1, appFct3)), , drop = FALSE]
  if (!identical(interval, "none")){
    SImat <- SImat[, -4, drop = FALSE]
    cNames <- c("Estimate", "Lower", "Upper")
  } else {
    cNames <- c("Estimate", "Std. Error", "t-value", "p-value")
  }
  colnames(SImat) <- cNames
  ciLabel <- switch(interval, delta = "Delta method", tfls = "To and from log scale", fls = "From log scale", fieller = "Fieller")
  eval(parse(text="drc:::resPrint(SImat, 'Estimated ratios of effect doses', interval, ciLabel, display = display)"))  
}


createsifctm <- function (edfct, logBase = NULL, fls = FALSE, indexMat, lenCoef, intgrid, intweights, rfinterval){
  if (is.null(edfct)) {
    stop("SI values cannot be calculated")
  } else {
    if (!fls) {
      if (is.null(logBase)) {
        "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, intgrid, intweights, rfinterval, ...) {
          ED1 <- edfct(parm1, pair[1], reference = reference, type = type, intgrid, intweights, rfinterval, ...)
          ED1v <- ED1[[1]]
          ED1d <- rep(0, lenCoef)
          ED1d[indexMat[, jInd]] <- ED1[[2]]
          ED2 <- edfct(parm2, pair[2], reference = reference, type = type, intgrid, intweights, rfinterval, ...)
          ED2v <- ED2[[1]]
          ED2d <- rep(0, lenCoef)
          ED2d[indexMat[, kInd]] <- ED2[[2]]
          SIpair <- ED1v/ED2v
          SIder <- (ED1d - SIpair * ED2d)/ED2v
          return(list(val = SIpair, der = SIder, der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
        }
      } else {
        "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, intgrid, intweights, rfinterval, ...) {
          ED1 <- edfct(parm1, pair[1], reference = reference, type = type, intgrid, intweights, rfinterval, ...)
          ED1v <- ED1[[1]]
          ED1d <- rep(0, lenCoef)
          ED1d[indexMat[, jInd]] <- ED1[[2]]
          ED2 <- edfct(parm2, pair[2], reference = reference, type = type, intgrid, intweights, rfinterval, ...)
          ED2v <- ED2[[1]]
          ED2d <- rep(0, lenCoef)
          ED2d[indexMat[, kInd]] <- ED2[[2]]
          SIpair <- logBase^(ED1v - ED2v)
          SIder <- SIpair * log(logBase) * (ED1d - ED2d)
          return(list(val = SIpair, der = SIder, der1 = (log(logBase) * logBase^ED1v) * ED1d, der2 = (log(logBase) * logBase^ED2v) * ED2d, valnum = logBase^ED1v, valden = logBase^ED2v))
        }
      }
    } else {
      "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, intgrid, intweights, rfinterval, ...) {
        ED1 <- edfct(parm1, pair[1], reference = reference, type = type, intgrid, intweights, rfinterval, ...)
        ED1v <- ED1[[1]]
        ED1d <- rep(0, lenCoef)
        ED1d[indexMat[, jInd]] <- ED1[[2]]
        ED2 <- edfct(parm2, pair[2], reference = reference, type = type, intgrid, intweights, rfinterval, ...)
        ED2v <- ED2[[1]]
        ED2d <- rep(0, lenCoef)
        ED2d[indexMat[, kInd]] <- ED2[[2]]
        SIpair <- ED1v - ED2v
        SIder <- ED1d - ED2d
        return(list(val = SIpair, der = SIder, der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
      }
    }
    return(sifct)
  }
}



siInnerm <- function (indPair, pVec, compMatch, object, indexMat, parmMat, varMat, level, reference, type, sifct, interval, degfree, logBase, intgrid, intweights, rfinterval){
  jInd <- indPair[1]
  kInd <- indPair[2]
  parmInd1 <- indexMat[, jInd]
  parmInd2 <- indexMat[, kInd]
  parmChosen1 <- parmMat[, jInd]
  parmChosen2 <- parmMat[, kInd]
  SIeval <- sifct(parmChosen1, parmChosen2, pVec, jInd, kInd, reference, type, intgrid, intweights, rfinterval)
  SIval <- SIeval$val
  dSIval <- SIeval$der
  oriMatRow <- c(SIval, sqrt(t(dSIval) %*% varMat %*% dSIval))
  siMatRow <- matrix(NA, 1, 4)
  siMatRow[1, 1] <- SIval
  if (identical(object$type, "continuous")){
    qFct <- function(x) {
      qt(x, degfree)
    }
    pFct <- function(x) {
      pt(x, degfree)
    }
  } else {
    qFct <- qnorm
    pFct <- pnorm
  }
  if (identical(interval, "none")){
    siMatRow[2] <- oriMatRow[2]
    tempStat <- (siMatRow[1] - 1)/siMatRow[2]
    siMatRow[3] <- tempStat
    siMatRow[4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))
  }
  if ((identical(interval, "delta")) || (identical(interval, "fls"))){
    stErr <- oriMatRow[2]
    tquan <- qFct(1 - (1 - level)/2)
    siMatRow[2] <- siMatRow[1] - tquan * stErr
    siMatRow[3] <- siMatRow[1] + tquan * stErr
    ciLabel <- "Delta method"
  }
  if (identical(interval, "tfls")){
    lsVal <- log(oriMatRow[1])
    lsdVal <- oriMatRow[2]/oriMatRow[1]
    tquan <- qFct(1 - (1 - level)/2)
    siMatRow[2] <- exp(lsVal - tquan * lsdVal)
    siMatRow[3] <- exp(lsVal + tquan * lsdVal)
    ciLabel <- "To and from log scale"
  }
  if ((!is.null(logBase)) && (identical(interval, "fls"))){
    siMatRow[1] <- logBase^(siMatRow[1])
    siMatRow[2] <- logBase^(siMatRow[2])
    siMatRow[3] <- logBase^(siMatRow[3])
    ciLabel <- "From log scale"
  }
  if (identical(interval, "fieller")){
    vcMat <- matrix(NA, 2, 2)
    vcMat[1, 1] <- SIeval$der1 %*% varMat[parmInd1, parmInd1] %*% SIeval$der1
    vcMat[2, 2] <- SIeval$der2 %*% varMat[parmInd2, parmInd2] %*% SIeval$der2
    vcMat[1, 2] <- SIeval$der1 %*% varMat[parmInd1, parmInd2] %*% SIeval$der2
    vcMat[2, 1] <- vcMat[1, 2]
    muVec <- c(SIeval$valnum, SIeval$valden)
    siMatRow[2:3] <- eval(parse(text="drc:::fieller(muVec, degfree, vcMat, level = level)"))
    ciLabel <- "Fieller"
  }
  siMatRow
}

