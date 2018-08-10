\name{summary-methods}
\docType{methods}
\alias{summary,IdtMclust-method}
\title{ Methods for Function summary in Package `MAINT.Data'}

\description{summary methods for classes defined in Package `MAINT.Data'.}

\usage{
 \S4method{summary}{IdtMclust}(object, parameters = FALSE, classification = FALSE, \dots)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{IdtMclust}} representing the results of fitting Gaussian mixtures to interval data objects}
  \item{parameters}{A boolean flag indicating if the parameter estimates of the optimal mixture should be displayed}
  \item{classification}{A boolean flag indicating if the crisp classification resulting from the optimal mixture should be displayed}
  \item{\dots}{Other named arguments.}
}

\keyword{methods}
\keyword{interval data}



