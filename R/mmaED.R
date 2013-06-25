mmaED <-
function(object, ..., respLev, ic = "AIC", interval = c("none", "buckland", "kang"), level = 0.9, bmd = c("none", "additional", "extra"), 
background = 0.05, dmList = NULL){
  lllist <- list(object, ...)
  ismedrc <- sapply(lllist, function(x) inherits(x, "medrc") | inherits(x, "glsdrc"))
  mllist <- lllist[ismedrc]  
  msl <- sapply(mllist, function(x) x$mselect)
  Call <- match.call()
  Call$respLev <- NULL
  Call$ic <- NULL
  Call$interval <- NULL
  Call$level <- NULL
  Call$bmd <- NULL
  Call$background <- NULL
  Call$dmList <- NULL
  colnames(msl) <- as.character(Call[-1L])[ismedrc]
  aic <- t(msl)[,ic]
  d <- aic-min(aic)
  wi <- exp(-0.5*d)/sum(exp(-0.5*d))
  names(wi) <- colnames(msl)

  if (bmd[1] == "extra") respLev <- respLev * (1-background)  
  if (bmd[1] == "none"){
    respLevMat <- matrix(respLev, nrow=length(lllist), ncol=length(respLev), byrow=TRUE)    
  } else {
    lenRL <- length(respLev)
    respLevMat <- matrix(sapply(1:lenRL, function(i) sapply(1:length(lllist),function(x){
      objectFit <- lllist[[x]]$fit
      cVal <- try(fixef(objectFit)[["c"]], silent = TRUE)
      if (inherits(cVal, "try-error")) {cVal <- 0}  # setting c if missing (not estimated)
      dVal <- fixef(objectFit)[["d"]]
      if (cVal > dVal) {tempVal <- cVal; cVal <- dVal; dVal <- tempVal}
      varcorr <- VarCorr(objectFit)
      100 * (qnorm(1 - background) - qnorm(1-(background+respLev[i]/100))) * as.numeric(varcorr[attr(varcorr, "dimnames")[[1]] %in% "Residual", 2]) / (dVal - cVal)
    })), ncol=lenRL)
  }
  ## Calculating ED values
  if (is.null(dmList) || interval[1] == "kang"){
    edEst <- t(sapply(1:length(mllist), function(x) ED(mllist[[x]], respLevMat[x,], display = FALSE)[,1]))  
    if (ncol(respLevMat) == 1 & ncol(edEst) == 1) edEst <- t(edEst)
  } else {
    edESSE <- lapply(1:ncol(respLevMat), function(i) sapply(1:length(dmList), function(x){deltaMethod(mllist[[x]]$fit, dmList[[x]](respLevMat[x,i]/100), vcov(mllist[[x]]$fit))}))
    edEst <- matrix(sapply(edESSE, function(x) as.numeric(x[1, , drop = FALSE])), nrow=length(dmList), ncol=length(respLev))
    edSe <- matrix(sapply(edESSE, function(x) as.numeric(x[2, , drop = FALSE])), nrow=length(dmList), ncol=length(respLev))
  }
  rownames(edEst) <- names(wi)
  edVec <- apply(edEst * wi, 2, sum)

  if (interval[1] == "none"){
    retMat <- cbind("Estimate" = edVec)
  }  
  if (interval[1] == "buckland"){
    if (is.null(dmList)){
      edSe <- t(sapply(1:length(mllist), function(x) ED(mllist[[x]], respLevMat[x,], display = FALSE)[,2]))
      if (ncol(respLevMat) == 1) edSe <- t(edSe)     
    }
    seVec <- apply(sqrt(edSe^2 + (t(t(edEst) - apply(edEst, 2, mean)))^2) * wi, 2, sum)
    quantVal <- qnorm(1 - (1 - level)/2) * seVec
    retMat <- as.matrix(cbind(edVec, seVec, edVec - quantVal, edVec + quantVal))
    colnames(retMat) <- c("Estimate", "SE", "Lower", "Upper")
  }

  if (interval[1] == "kang"){
    edCll <- t(sapply(1:length(mllist), function(x) ED(mllist[[x]], respLevMat[x,], interval="delta", level = level, display = FALSE)[,3]))
    edClu <- t(sapply(1:length(mllist), function(x) ED(mllist[[x]], respLevMat[x,], interval="delta", level = level, display = FALSE)[,4]))
    if (ncol(respLevMat) == 1){
      edCll <- t(edCll)
      edClu <- t(edClu)
    }
    retMat <- as.matrix(cbind(edVec, apply(edCll * wi, 2, sum), apply(edClu * wi, 2, sum)))
    colnames(retMat) <- c("Estimate", "Lower", "Upper")
  }  
  out <- list()
  out$effectiveRespLev <- if (bmd[1] == "none") respLev else respLevMat
  out$retMat <- retMat
  out$weights <- wi
  out$EDi <- edEst
  if (interval[1] == "buckland") {out$SEi <- edSe}  
  return(out)
}


