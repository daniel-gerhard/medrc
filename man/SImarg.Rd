\name{SImarg}

\alias{SImarg}
\title{Comparing marginal selectivity indices across curves}

\description{
  'SI' compares selectivity indices for arbitrary dosage across curves.
  The selectivity is the ratio between effective dosages from different curves. 
}

\usage{
  SImarg(object, percVec, compMatch = NULL, reverse = FALSE, 
  interval = c("none", "delta", "fieller", "fls"), 
  level = ifelse(!(interval == "none"), 0.95, NULL),
  reference = c("control", "upper"),  
  type = c("relative", "absolute"),
  nGQ=5, rfinterval=c(-1000, 1000),
  display = TRUE, logBase=NULL, ...)
}

\arguments{
  \item{object}{an object of class 'medrc'.}
  \item{percVec}{a numeric vector of dosage values.}
  \item{compMatch}{an optional character vector of names of assays to be compared. If not specified all comparisons are supplied.}
  \item{reverse}{logical. If TRUE the order of comparison of two curves is reversed.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fieller" for confidence intervals based on Fieller's theorem (with help from the delta method).
  Use "fls" for from-logarithm-scale-based confidence intervals (in case the parameter in the model fit is 
  log(ED50) as is the case for the \code{\link{logistic}} or \code{\link{llogistic2}} models); currently
  the argumetn \code{logBase} then also needs to be specified.}
  \item{level}{numeric. The level for the confidence intervals. Default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}  
  \item{type}{character string specifying whether absolute or relative response levels are supplied.}
  \item{nGQ}{integer. Specifies the number nof nodes for Gauss-Hermite quadrature.}
  \item{rfinterval}{numeric vector. Interval for root finding (uniroot) to search for ED values. Default is [-1000, 1000].}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{...}{additional arguments to the function doing the calculations.} 
}

\value{
  A matrix with columns containing the estimates, estimated standard errors, 
  t-statistics for testing indices equal to 1 and the corresponding p-values.
}

\references{
  Finney, D. J. (1978) \emph{Statistical method in Biological Assay}, London: Charles Griffin House, 
  3rd edition (pp. 80--82).

  Kotz, S. and Johnson, N. L. (1983) \emph{Encyclopedia of Statistical Sciences Volume 3}, 
  New York: Wiley \& Sons (pp. 86--87).
  
  Ritz, C. and Cedergreen, N. and Jensen, J. E. and Streibig, J. C. (2006) 
  Relative potency in nonsimilar dose-response curves, \emph{Weed Science}, \bold{54}, 407--412.
}

\author{Daniel Gerhard, Christian Ritz}


\keyword{models}
\keyword{nonlinear}
