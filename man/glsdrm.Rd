\name{glsdrm}
\alias{glsdrm}

\title{Dose-response curve estimation by generalized least squares}
\description{
Implements drc nonlinear functions into the nlme framework for nonlinear GLS dose-response modeling.
}
\usage{
glsdrm(form, curveid=NULL, data, fct, correlation = NULL, 
       weights = NULL, control = NULL, start = NULL)
}

\arguments{
  \item{form}{Formula describing the dose-response relationship}
  \item{curveid}{Formula with parameter names on the left hand side (divided by +) and a column name in data, denoting a factor, to estimate separate parameters per factor-level. If NULL only fixed effects for a single curve will be estimated.}
  \item{data}{a data.frame object}
  \item{fct}{a function of the drc package}
  \item{correlation}{additional corClasses object}
  \item{weights}{additional varClasses object}
  \item{control}{list with nlme control arguments}
  \item{start}{optional list with initial values for the fixed components. If NULL the initial values will be found automatically.}
}

\value{An object of class glsdrc}

\author{Daniel Gerhard}

\seealso{
  \code{\link{drm}}, \code{\link{nlme}}
}

\keyword{ models }
\keyword{ nonlinear }

