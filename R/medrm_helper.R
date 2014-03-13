makehelpfunction <- function(fct){
  eval(parse(text="drcfunction <- function(){
    if (is.null(fct$fixed)) fct$fixed <- rep(NA, length(fct$names))
    parmVec <- fct$fixed
    notFixed <- is.na(parmVec)
    numParm <- length(parmVec)
    .value <- fct$fct(dose, eval(parse(text=paste('cbind(', paste(fct$names, collapse=', '),')', sep = ''))))
    .actualArgs <- as.list(match.call()[fct$names])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
      .grad <- fct$deriv1(dose, eval(parse(text=paste('cbind(', paste(fct$names, collapse=', '),')', sep = ''))))
      dimnames(.grad) <- list(NULL, .actualArgs)
      attr(.value, 'gradient') <- .grad
    }
    .value
  }"))
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
    fullp <- coefficients(drm(form, data=mf[,1:2], fct=fct))
    pmat <- t(sapply(spl, function(x){
      trycoef <- try(coefficients(drm(form, data=x, fct=fct)), silent=TRUE)  
      if (class(trycoef)[1] == "try-error") return(fullp) else return(trycoef)
    }))     
    plist <- lapply(1:length(fid), function(i){
      if (fid[i]) fullp[i] else pmat[,i]   
    }) 
    return(unlist(plist))
  }
}


findrandomstartvals <- function(form, data, fct, random){
  if (class(random)[1] == "list"){
    rid <- names(random)
    rpars <- lapply(random, function(x){
      spr <- strsplit(deparse(x[[2]]), "+")[[1]]
      return(spr[!(spr %in% c(" ", "+"))])
    })
  }
  if (class(random)[1] == "formula"){
    rid <- as.character(random[[3]][[3]])
    rid <- rid[rid != "/"]
    spr <- strsplit(deparse(random[[2]]), "+")[[1]]
    rpars <- list()
    for (i in 1:length(rid)){
      rpars[[i]] <- spr[!(spr %in% c(" ", "+"))]
    }
  }
  rform <- as.formula(paste(as.character(form)[2], as.character(form)[1], as.character(form)[3], "+", paste(rid, collapse="+")))
  mf <- model.frame(rform, data)
  # 2nd hierarchy compared with fullp instead of first hierarchical level
  rstart <- lapply(1:length(rid), function(ri){
    sf1 <- apply(mf[,rid[1:ri],drop=FALSE], 1, function(x) paste(x, collapse="/"))
    spl <- split(mf, sf1)
    if (ri == 1){
      fullpp <- coefficients(drm(form, data=mf[,1:2], fct=fct))
      fullp <- matrix(rep(fullpp, each=length(spl)), ncol=length(fullpp))
    } else {
      sf2 <- apply(mf[,rid[1:(ri-1)],drop=FALSE], 1, function(x) paste(x, collapse="/"))
      spl2 <- split(mf, sf2)
      fullpp <- coefficients(drm(form, data=mf[,1:2], fct=fct))
      fullp <- t(sapply(spl2, function(x){
        trycoef <- try(coefficients(drm(form, data=x, fct=fct)), silent=TRUE)  
        if (class(trycoef)[1] == "try-error") return(fullp) else return(trycoef)
      })) 
      fullp <- matrix(unlist(lapply(1:length(spl2), function(c2){
        nrep <- length(unique(apply(spl2[[c2]][,rid[1:ri],drop=FALSE], 1, function(x) paste(x, collapse="/"))))
        matrix(rep(fullp[c2,], nrep), ncol=nrep)
      })), ncol=length(fullpp), byrow=TRUE)
    }
    
    pmat <- t(sapply(spl, function(x){
      trycoef <- try(coefficients(drm(form, data=x, fct=fct)), silent=TRUE)  
      if (class(trycoef)[1] == "try-error") return(fullpp) else return(trycoef)
    }))     
    
    rmat <- pmat-fullp
    colnames(rmat) <- fct$names
    return(rmat)
  })
  names(rstart) <- rid
  
  for (i in 1:length(rid)){
    rownames(rstart[[i]]) <- unique(apply(mf[,rid[1:i],drop=FALSE], 1, function(x) paste(x, collapse="/")))
    rstart[[i]] <- rstart[[i]][,rpars[[i]], drop=FALSE]
  }    
  return(rstart)  
}
