\name{mmaED}
\alias{mmaED}
\title{ED estimation for averaged medrc models}
\description{Estimates ED for averaged fixed effects of several medrc objects}
\usage{
mmaED(object, ..., respLev, ic = "AIC", interval = c("none", "buckland", "kang"), 
     level = 0.9, bmd = c("none", "additional", "extra"), 
     background = 0.05, dmList = NULL)
}
\arguments{
  \item{object}{an medrc object}
  \item{\dots}{further medrc objects}
  \item{respLev}{a numeric vector containing the response levels}
  \item{ic}{the information criterion used for model averaging}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{level}{confidence level, default at 0.9}
  \item{bmd}{benchmark dose estimation (smallest dose resulting in a
    probability of an abnormal response)}
  \item{background}{probability of an abnormal response}
  \item{dmList}{list of functions, estimating bmd on a transformed scale}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
