\name{cor--methods}
\docType{methods}
\alias{cor}
\alias{cor-methods}
\alias{cor,IdtNDE-method}
\alias{cor,IdtSNDE-method}
\alias{cor,IdtNandSNDE-method}
\alias{cor,IdtMxNDE-method}
\alias{cor,IdtMxSNDE-method}

\title{Methods for function cor in Package `MAINT.Data'}
\description{S4 methods for function cor. These methods extract estimates of correlation matrices for the models fitted to Interval Data.}

\keyword{methods}
\keyword{cor}
\keyword{Interval Data}

\usage{

\S4method{cor}{IdtNDE}(x)
\S4method{cor}{IdtSNDE}(x)
\S4method{cor}{IdtNandSNDE}(x)
\S4method{cor}{IdtMxNDE}(x)
\S4method{cor}{IdtMxSNDE}(x)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}

}

\value{For the \code{\linkS4class{IdtNDE}}, \code{\linkS4class{IdtSNDE}} and \code{\linkS4class{IdtNandSNDE}} methods or \code{\linkS4class{IdtMxNDE}}, \code{\linkS4class{IdtMxSNDE}} methods with slot \dQuote{Hmcdt} equal to TRUE: a matrix with the estimated correlations. \cr For the \code{\linkS4class{IdtMxNDE}}, and \code{\linkS4class{IdtMxSNDE}} methods with slot \dQuote{Hmcdt} equal to FALSE: a three-dimensional array with a matrix with the estimated correlations for each group at each level of the third dimension.}

\seealso{
  \code{\link{var}} 
}

