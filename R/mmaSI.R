mmaSI <- function(object, ..., percVec, ic = c("AIC", "BIC"), interval = c("none", "buckland", "kang"), level = 0.95, compMatch = NULL, reverse = FALSE, marginal=FALSE, nGQ=5, rfinterval=c(-1000, 1000)){
  interval <- match.arg(interval)
  ic <- match.arg(ic)
  lllist <- list(object, ...)
  ismedrc <- sapply(lllist, function(x) inherits(x, "medrc") | inherits(x, "glsdrc"))
  mllist <- lllist[ismedrc]  
  Call <- match.call()
  Call$percVec <- NULL
  Call$ic <- NULL
  Call$interval <- NULL
  Call$level <- NULL
  Call$compMatch <- NULL
  Call$reverse <- NULL
  Call$marginal <- NULL
  Call$nGQ <- NULL
  Call$rfinterval <- NULL
  mnames <- as.character(Call[-1L])[ismedrc]
  # number of curves
  if (any(ncol(mllist[[1]]$parmMat) != sapply(mllist, function(x) ncol(x$parmMat)))) stop("Number of curves are not the same for all models!")
  # model weights
  if (ic == "AIC") inc <- sapply(mllist, AIC)
  if (ic == "BIC") inc <- sapply(mllist, BIC)
  d <- inc-min(inc)
  wi <- exp(-0.5*d)/sum(exp(-0.5*d))
  names(wi) <- mnames
  
  # calculate ED per model
  if (marginal == TRUE){
    silist <- lapply(mllist, function(x) SImarg(x, percVec=percVec, display=FALSE, compMatch=compMatch, reverse=reverse, nGQ=nGQ, rfinterval=rfinterval))
  } else {
    silist <- lapply(mllist, function(x) SI(x, percVec=percVec, display=FALSE, compMatch=compMatch, reverse=reverse))
  }
  sim <- rbind(sapply(silist, function(x) x[,1]))
  colnames(sim) <- mnames
  
  # model-averaged ED
  sima <- apply(t(t(sim) * wi), 1, sum)
  
  if (interval[1] == "none"){
    retMat <- cbind("Estimate" = sima)
  }  
  
  if (interval[1] == "buckland"){
    sise <- rbind(sapply(silist, function(x) x[,2]))
    sema <- apply(sqrt(sise^2 + (sim - apply(sim, 1, mean))^2), 1, function(x) sum(x*wi))
    quant <- qnorm(1 - (1 - level)/2) * sema
    retMat <- as.matrix(cbind(sima, sema, sima - quant, sima + quant))
    colnames(retMat) <- c("Estimate", "SE", "Lower", "Upper")
  }
  
  if (interval[1] == "kang"){
    if (marginal == TRUE){
      siclist <- lapply(mllist, function(x) SImarg(x, percVec=percVec, interval="delta", display=FALSE, nGQ=nGQ, rfinterval=rfinterval))
    } else {
      siclist <- lapply(mllist, function(x) SI(x, percVec=percVec, interval="delta", display=FALSE))
    }
    ll <- rbind(sapply(siclist, function(x) x[,2]))
    ul <- rbind(sapply(siclist, function(x) x[,3]))
    retMat <- as.matrix(cbind(sima, apply(ll, 1, function(x) sum(x*wi)), apply(ul, 1, function(x) sum(x*wi))))
    colnames(retMat) <- c("Estimate", "Lower", "Upper")
  }  
  
  rownames(retMat) <- rownames(silist[[1]])
  
  out <- list()
  out$retMat <- retMat
  out$weights <- wi
  out$SIi <- sim
  if (interval[1] == "buckland") {out$SEi <- sise}  
  return(out)
}
