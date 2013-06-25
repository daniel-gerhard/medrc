residuals.medrc <- function(object, ...){
  fct <- object$fct
  makehelpfunction(fct)
  resid(object$fit, type="n")
}

predict.medrc <- function(object, ..., newdata=NULL, level=NULL){
  fct <- object$fct
  makehelpfunction(fct)
  if (is.null(newdata)) newdata <- object$data
  if (is.null(level)) level <- object$fit$dims$Q
  predict(object$fit, newdata=newdata, level=level)
}


ranef.medrc <- function(object, ...){
  ranef(object$fit, ...)
}


vcov.medrc <- function(object, ...){
  #est <- fixef(object$fit)
  #stdFixed <- sqrt(diag(as.matrix(object$fit$varFix)))
  #std <- sqrt(object$fit$dims$N/(object$fit$dims$N - length(stdFixed))) * stdFixed
  #cr <- array(t(object$fit$varFix/stdFixed)/stdFixed, dim(object$fit$varFix), list(names(est), names(est)))
  #return(cor2cov(round(cr,10), std))
  object$vc
}

summary.medrc <- function(object, ...){
  summary(object$fit)
}


df.residual.medrc <- function(object, ...){
  # need to define a better residual.df !!!!!
  object$fit$dims$N - length(coefficients(object)) - nrow(VarCorr(object$fit))
}


AIC.medrc <- function(object, ..., k = 2){
  objlist <- list(object, ...)
  if (length(objlist) > 1){
    ismedrc <- sapply(objlist, function(x) inherits(x, "medrc"))
    medrclist <- objlist[ismedrc]  
    Call <- match.call()
    Call$k <- NULL
    names(medrclist) <- as.character(Call[-1L])[ismedrc]
    nlmelist <- lapply(medrclist, function(x) x$fit)
    ftext <- paste("AIC(", paste(paste("nlmelist[['", names(nlmelist), "']]", sep=""), collapse=","), ", k=", k, ")", sep="")
    AICtab <- eval(parse(text=ftext))
    rownames(AICtab) <- names(nlmelist)
  } else {
    AICtab <- AIC(object$fit, k=k)
  }
  return(AICtab)
}

BIC.medrc <- function(object, ...){
  objlist <- list(object, ...)
  if (length(objlist) > 1){
    ismedrc <- sapply(objlist, function(x) inherits(x, "medrc"))
    medrclist <- objlist[ismedrc]  
    Call <- match.call()
    Call$k <- NULL
    names(medrclist) <- as.character(Call[-1L])[ismedrc]
    nlmelist <- lapply(medrclist, function(x) x$fit)
    ftext <- paste("BIC(", paste(paste("nlmelist[['", names(nlmelist), "']]", sep=""), collapse=","), ")", sep="")
    BICtab <- eval(parse(text=ftext))
    rownames(BICtab) <- names(nlmelist)
  } else {
    BICtab <- BIC(object$fit)
  }
  return(BICtab)
}

logLik.medrc <- function(object, REML = FALSE, ...){
  logLik(object$fit, REML=REML, ...) 
}



