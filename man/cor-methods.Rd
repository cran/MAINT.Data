\name{cor--methods}
\docType{methods}
\alias{cor}
\alias{cor-methods}
\alias{cor,IdtNDE-method}
\alias{cor,IdtSNDE-method}
\alias{cor,IdtNandSNDE-method}

\title{Methods for function cor in Package `MAINT.Data'}
\description{S4 methods for function cor. These methods extract estimates of correlation matrices for the models fitted to Interval Data.}

\keyword{methods}
\keyword{cor}
\keyword{Interval Data}

\usage{

\S4method{cor}{IdtNDE}(x)
\S4method{cor}{IdtSNDE}(x)
\S4method{cor}{IdtNandSNDE}(x)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}

}

\value{For the \dQuote{IdtNDE}, \dQuote{IdtSNDE} and \dQuote{IdtSNDE} methods or \dQuote{IdtMxNDE}, \dQuote{IdtMxSNDE} methods with slot \dQuote{Hmcdt} equal to TRUE: a matrix of the estimated correlations. For the \dQuote{IdtMxNDE}, and \dQuote{IdtMxSNDE} methods with slot \dQuote{Hmcdt} equal to FALSE: a three-dimensional array with a matrix of the estimated correlations for each group at each level of the third dimension.}

\seealso{
  \code{\link{var}} 
}

