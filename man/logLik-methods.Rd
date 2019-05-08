\name{LogLikelihood--methods}
\docType{methods}
\alias{logLik}
\alias{logLik-methods}
\alias{logLik,IdtE-method}

\title{Methods for function logLik in Package `MAINT.Data'}
\description{S4 methods for function logLik. These methods extract the logarithm of the likelihood function for the models fitted to Interval Data.}

\keyword{methods}
\keyword{logLik}
\keyword{Interval Data}

\usage{

\S4method{logLik}{IdtE}(x)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}
}

\value{The information criterion value for the model fitted in object x.}

\seealso{
  \linkS4class{IdtE}, \code{\link{BIC}}, \code{\link{AIC}}  
}

