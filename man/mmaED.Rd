\name{mmaED}
\alias{mmaED}
\title{ED estimation for averaged medrc and glsdrc models}
\description{Estimates ED for averaged fixed effects of several medrc or glsdrc objects}
\usage{
mmaED(object, ..., respLev, ic = "AIC", interval = c("none", "buckland", "kang"), 
     level = 0.95)
}
\arguments{
  \item{object}{an medrc object}
  \item{\dots}{further medrc objects}
  \item{respLev}{a numeric vector containing the response levels}
  \item{ic}{the information criterion used for model averaging}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{level}{confidence level, default at 0.95}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
