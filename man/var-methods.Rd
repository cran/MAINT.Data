\name{var--methods}
\docType{methods}
\alias{var}
\alias{var-methods}
\alias{var,IdtNDE-method}
\alias{var,IdtSNDE-method}
\alias{var,IdtNandSNDE-method}

\title{Methods for function var in Package `MAINT.Data'}
\description{S4 methods for function var. These methods extract estimates of variance-covariance matrices for the models fitted to Interval Data.}

\keyword{methods}
\keyword{var}
\keyword{Interval Data}

\usage{

\S4method{var}{IdtNDE}(x)
\S4method{var}{IdtSNDE}(x)
\S4method{var}{IdtNandSNDE}(x)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}

}

\value{For the \dQuote{IdtNDE}, \dQuote{IdtSNDE} and \dQuote{IdtSNDE} methods or \dQuote{IdtMxNDE}, \dQuote{IdtMxSNDE} methods with slot \dQuote{Hmcdt} equal to TRUE: a matrix of the estimated covariances. For the \dQuote{IdtMxNDE}, and \dQuote{IdtMxSNDE} methods with slot \dQuote{Hmcdt} equal to FALSE: a three-dimensional array with a matrix of the estimated covariances for each group at each level of the third dimension.}


\seealso{
  \code{\link{cor}} 
}

