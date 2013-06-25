makehelpfunction <- function(fct){
  drcfunction <- function(){
    if (is.null(fct$fixed)) fct$fixed <- rep(NA, length(fct$names))
    parmVec <- fct$fixed
    notFixed <- is.na(parmVec)
    numParm <- length(parmVec)
    .value <- fct$fct(dose, eval(parse(text=paste("cbind(", paste(fct$names, collapse=", "),")", sep = ""))))
    .actualArgs <- as.list(match.call()[fct$names])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
      .grad <- fct$deriv1(dose, eval(parse(text=paste("cbind(", paste(fct$names, collapse=", "),")", sep = ""))))
      dimnames(.grad) <- list(NULL, .actualArgs)
      attr(.value, "gradient") <- .grad
    }
    .value
  }
  aargs <- paste("dose=,", paste(fct$names, collapse="=,"), "=", sep="")
  formals(drcfunction) <- eval(parse(text = paste("alist(", aargs,")", sep = "")))
  attr(drcfunction, "pnames") <- fct$names
  attr(drcfunction, "class") <- "selfStart"
  attr(drcfunction,"initial") <- function(mCall, data, LHS){
    xy <- sortedXyData(mCall[["dose"]], LHS, data)
    val <- fct$ssfct(xy)
    names(val) <- fct$names
    val
  }
  assign("drcfunction", drcfunction, envir=.GlobalEnv) # find any way around assigning to global envir??? 
}



findfixedstartvals <- function(form, data, cid, fct, fid, mform){
  if (is.na(cid)){
    mf <- model.frame(form, data)
    return(fct$ssfct(mf))
  } else {
    cform <- as.formula(paste(as.character(form)[2], as.character(form)[1], as.character(form)[3], "+", cid))
    mf <- model.frame(cform, data)
    spl <- split(mf, mf[,3])
    fullp <- coefficients(eval(parse(text=paste("nls(", mform, ", data=mf[,1:2])", sep=""))))
    pmat <- t(sapply(spl, function(x){
      trycoef <- try(coefficients(eval(parse(text=paste("nls(", mform, ", data=x)", sep="")))), silent=TRUE)  
      if (class(trycoef)[1] == "try-error") return(fullp) else return(trycoef)
    }))   
    plist <- lapply(1:length(fid), function(i){
      if (fid[i]) fullp[i] else pmat[,i]   
    }) 
    return(unlist(plist))
  }
}
