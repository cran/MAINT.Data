\name{coef--methods}
\docType{methods}
\alias{coef}
\alias{coef-methods}
\alias{coef,IdtNDE-method}
\alias{coef,IdtSNDE-method}
\alias{coef,IdtNandSNDE-method}

\title{Methods for function coef in Package `MAINT.Data'}
\description{S4 methods for function coef. As in the generic coef S3 `stats' method, these methods extract parameter estimates for the models fitted to Interval Data.}

\keyword{methods}
\keyword{coef}
\keyword{Interval Data}

\usage{

\S4method{coef}{IdtNDE}(object, selmodel=BestModel(object), \dots)
\S4method{coef}{IdtSNDE}(object, selmodel=BestModel(object), ParType=c("Centr", "Direct", "All"), \dots)
\S4method{coef}{IdtNandSNDE}(object, selmodel=BestModel(object),  ParType=c("Centr", "Direct", "All"), \dots)

}

\arguments{
  \item{object}{ An object representing a model fitted to interval data.}
  \item{selmodel}{Selected model from a list of candidate models saved in object.}
  \item{ParType}{Parameterization of the Skew-Normal distribution. Alternatives are \dQuote{Centr} for centred parameters,
\dQuote{Direct} for direct parameters and \dQuote{All}, for both types of parameters. See Arellano-Valle and Azzalini (2008) for details.}
  \item{\dots}{Additional arguments for method functions.}
}

\value{ A vector of parameter estimates.}

\references{
Arellano-Valle, R. B. and Azzalini, A. (2008): "The centred parametrization for the multivariate skew-normal distribution". Journal of Multivariate Analysis, Volume 99, Issue 7, 1362-1382.
}

\seealso{
  \code{\link{stdEr}}, \code{\link{vcov}} 
}

