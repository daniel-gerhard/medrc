\name{mmaBMD}
\alias{mmaBMD}
\title{BMD estimation for averaged medrc or glsdrc models}
\description{Estimates the benchmark dose for averaged fixed effects of several medrc or glsdrc objects}
\usage{
mmaBMD(object, ..., respLev, ic = "AIC", interval = c("none", "buckland", "kang"), 
     level = 0.95, bmd = c("additional", "extra"), 
     background = 0.05)
}
\arguments{
  \item{object}{an medrc object}
  \item{\dots}{further medrc objects}
  \item{respLev}{a numeric vector containing the benchmark response levels}
  \item{ic}{the information criterion used for model averaging}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{level}{confidence level, default at 0.95}
  \item{bmd}{benchmark dose estimation (smallest dose resulting in a
    probability of an abnormal response)}
  \item{background}{probability of an abnormal response}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
