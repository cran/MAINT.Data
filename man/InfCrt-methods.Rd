\name{InformationCriteria--methods}
\docType{methods}
\alias{BIC}
\alias{AIC}
\alias{BIC-method}
\alias{AIC-method}
\alias{BIC,IdtE-method}
\alias{AIC,IdtE-method}

\title{Methods for functions BIC and AIC in Package `MAINT.Data'}
\description{S4 methods for functions BIC and AIC. These methods extract the values of the BIC and AIC information criteria for the models fitted to Interval Data.}

\keyword{methods}
\keyword{BIC}
\keyword{AIC}
\keyword{Interval Data}

\usage{

\S4method{BIC}{IdtE}(x)
\S4method{AIC}{IdtE}(x)

}

\arguments{
  \item{x}{ An object representing a model fitted to interval data.}
}

\value{The information criterion value for the model fitted in object x.}

\seealso{
  \linkS4class{IdtE}, \code{\link{logLik}} 
}

