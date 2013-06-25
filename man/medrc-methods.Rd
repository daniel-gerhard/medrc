\name{medrc-methods}
\alias{df.residual.medrc}
\alias{predict.medrc}
\alias{vcov.medrc}
\alias{summary.medrc}
\alias{residuals.medrc}
\alias{ranef.medrc}
\alias{AIC.medrc}
\alias{BIC.medrc}
\alias{logLik.medrc}

\title{medrc methods}
\description{Methods for objects of class medrc}
\usage{
\method{df.residual}{medrc}(object, ...)
\method{predict}{medrc}(object, ..., newdata = NULL, level = NULL)
\method{residuals}{medrc}(object, ...)
\method{vcov}{medrc}(object, ...)
\method{summary}{medrc}(object, ...)
\method{ranef}{medrc}(object, ...)
\method{AIC}{medrc}(object, ..., k = 2)
\method{BIC}{medrc}(object, ...)
\method{logLik}{medrc}(object, REML = FALSE, ...)
}

\arguments{
  \item{object}{An object of class medrc}
  \item{newdata}{an optional data frame to be used for obtaining the predictions. All variables used in the fixed and random effects models, as well as the grouping factors, must be present in the data frame. If missing, the fitted values are returned.}
  \item{level}{an optional integer vector giving the level(s) of grouping to be used in obtaining the predictions. Level values increase from outermost to innermost grouping, with level zero corresponding to the population predictions. Defaults to the highest or innermost level of grouping.}
  \item{k}{numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
  \item{REML}{an optional logical value. If TRUE the restricted log-likelihood is returned, else, if FALSE, the log-likelihood is returned. Defaults to FALSE.}
  \item{...}{further arguments}
}

\author{Daniel Gerhard}



\keyword{ methods }
