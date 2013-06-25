residuals.glsdrc <- function(object, ...){
  fct <- object$fct
  makehelpfunction(fct)
  resid(object$fit, type="n")
}

predict.glsdrc <- function(object, ..., newdata=NULL){
  fct <- object$fct
  makehelpfunction(fct)
  if (is.null(newdata)) newdata <- object$data
  predict(object$fit, newdata=newdata)
}

vcov.glsdrc <- function(object, ...){
  return(object$vc)
}

summary.glsdrc <- function(object, ...){
  summary(object$fit)
}

df.residual.glsdrc <- function(object, ...){
  object$fit$dims$N - object$fit$dims$p
}

AIC.glsdrc <- function(object, ..., k = 2){
  objlist <- list(object, ...)
  if (length(objlist) > 1){
    isglsdrc <- sapply(objlist, function(x) inherits(x, "glsdrc"))
    glsdrclist <- objlist[isglsdrc]  
    Call <- match.call()
    Call$k <- NULL
    names(glsdrclist) <- as.character(Call[-1L])[isglsdrc]
    nlmelist <- lapply(glsdrclist, function(x) x$fit)
    ftext <- paste("AIC(", paste(paste("nlmelist[['", names(nlmelist), "']]", sep=""), collapse=","), ", k=", k, ")", sep="")
    AICtab <- eval(parse(text=ftext))
    rownames(AICtab) <- names(nlmelist)
  } else {
    AICtab <- AIC(object$fit, k=k)
  }
  return(AICtab)
}

BIC.glsdrc <- function(object, ...){
  objlist <- list(object, ...)
  if (length(objlist) > 1){
    isglsdrc <- sapply(objlist, function(x) inherits(x, "glsdrc"))
    glsdrclist <- objlist[isglsdrc]  
    Call <- match.call()
    Call$k <- NULL
    names(glsdrclist) <- as.character(Call[-1L])[isglsdrc]
    nlmelist <- lapply(glsdrclist, function(x) x$fit)
    ftext <- paste("BIC(", paste(paste("nlmelist[['", names(nlmelist), "']]", sep=""), collapse=","), ")", sep="")
    BICtab <- eval(parse(text=ftext))
    rownames(BICtab) <- names(nlmelist)
  } else {
    BICtab <- BIC(object$fit)
  }
  return(BICtab)
}

logLik.glsdrc <- function(object, REML = FALSE, ...){
  logLik(object$fit, REML=REML, ...) 
}

