\name{EDcond}

\alias{EDcond}


\title{Estimating effective doses conditional on random effects}

\description{

  \code{ED} estimates effective doses (ECp/EDp/ICp) for given reponse levels.
}

\usage{
  EDcond(object, respLev, cranef=NULL, interval = c("none", "delta", "fls", "tfls"), 
         level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), 
         type = c("relative", "absolute"), lref, uref, bound = TRUE, 
         display = TRUE, logBase = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'medrc'.}
  \item{respLev}{a numeric vector containing the response levels.}
  \item{cranef}{a matrix containing random effects for conditional ED estimation. Column names need to be the same as the model parameters (b,c,d,e,...). If NULL, the random effects of the model object are used. Only works with a single variance component and single intercepts per non-linear model parameters.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fls" for from logarithm scale based confidence intervals (in case the parameter in the model is log(ED50) as for
  the \code{\link{llogistic2}}) models. The only alternative for model-robust fits is using inverse regression.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{lref}{numeric value specifying the lower limit to serve reference.}
  \item{uref}{numeric value specifying the upper limit to serve reference (eg. 100\%).}  
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}  
  \item{...}{see the details section below.}
}


\value{
  A matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz, Daniel Gerhard}

%\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\keyword{models}
\keyword{nonlinear}
