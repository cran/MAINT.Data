\name{stdEr}
\alias{stdEr}

\title{Standard deviations}
\description{
   Extract standard deviations from estimated models.
}
\usage{
  stdEr(x, \dots)
}

\arguments{
  \item{x}{an R object containing estimates of a statistical model}
  \item{\dots}{further arguments for methods}
}

\details{
  \code{stdEr} is a generic function. The default method returns the square root of the diagonal of
  the variance-covariance matrix.
}

\value{
  numeric, the estimated standard errors of the parameter model estimates.  
}

\seealso{\code{\link{vcov}}, \code{\link{summary}}.}

\keyword{methods}
