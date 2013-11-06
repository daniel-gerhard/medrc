\name{mmaSI}
\alias{mmaSI}
\title{Comparing model averaged selectivity indices across curves for medrc and glsdrc models}
\description{Compares selectivity indices for model averaged fixed effects of several medrc or glsdrc objects}
\usage{
mmaSI(object, ..., percVec, ic = c("AIC", "BIC"), interval = c("none", "buckland", "kang"), 
     level = 0.95, compMatch = NULL, reverse=FALSE, marginal=FALSE, nGQ=5, rfinterval=c(-1000, 1000))
}
\arguments{
  \item{object}{an medrc object}
  \item{\dots}{further medrc objects}
  \item{percVec}{a numeric vector of response levels}
  \item{ic}{the information criterion used for model averaging}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{level}{confidence level, default at 0.95}
  \item{compMatch}{an optional character vector of names of assays to be compared. If not specified all comparisons are supplied}
  \item{reverse}{logical. If TRUE the order of comparison of two curves is reversed}
  \item{marginal}{logical. If TRUE, marginalized SI estimates are obtained. See SImarg for further information.}
  \item{nGQ}{integer. Specifies the number nof nodes for Gauss-Hermite quadrature.}
  \item{rfinterval}{numeric vector. Interval for root finding (uniroot) to search for ED values. Default is [-1000, 1000].}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
