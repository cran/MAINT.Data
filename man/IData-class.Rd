\name{IData-class}
\Rdversion{1.1}
\docType{class}
\alias{IData-class}
\alias{head,IData-method}
\alias{print,IData-method}
\alias{summary,IData-method}
\alias{tail,IData-method}

\title{Class "IData"}
\description{
  A data-array of interval-valued data is an array where each of the NObs rows, corresponding to each entity under analysis, contains the observed intervals of the NIVar descriptive variables.
}

\section{Slots}{
  \describe{
    \item{\code{MidP}:}{A data-frame of the midpoints of the observed intervals}
    \item{\code{LogR}:}{A data-frame of the logarithms of the ranges of the observed intervals}
    \item{\code{ObsNames}:}{An optional vector of names assigned to the individual observations.}
    \item{\code{VarNames}:}{An optional vector of names to be assigned to the Interval-Valued Variables.}
    \item{\code{NObs}:}{Number of entities under analysis (cases)}
    \item{\code{NIVar}:}{Number of interval variables}
  }
}

\section{Methods}{
  \describe{
    \item{head}{\code{signature(x = "IData")}: head S4 method for the IData-class}
    \item{show}{\code{signature(object = "IData")}: show S4 method for the IData-class}
    \item{nrow}{\code{signature(x = "IData")}: returns the number of entities}
    \item{ncol}{\code{signature(x = "IData")}: returns the number of Iterval Variables }
    \item{tail}{\code{signature(x = "IData")}: tail S4 method for the IData-class }
    \item{mle}{\code{signature(x = "IData")}: Maximum likelihood estimation}
    \item{MANOVA}{\code{signature(x = "IData")}: MANOVA tests on the interval-valued data}
    \item{lda}{\code{signature(x = "IData")}: performs Linear Discriminant Analysis using the estimated model parameters}
    \item{qda}{\code{signature(x = "IData")}: performs Quadratic Discriminant Analysis using the estimated model parameters} }
}

\references{
Brito, P., Duarte Silva, A. P.  (2011): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics (forthcoming).\cr

Noirhomme-Fraiture, M., Brito, P. (2011): "Far Beyond the Classical Data Models: Symbolic Data Analysis". Statistical Analysis and Data Mining, Volume 4, Issue 2, 157-170.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\seealso{
\code{\link{IData}}, \code{\link{mle}}, \code{\link{MANOVA}}, \code{\link{lda}}, \code{\link{qda}}
}

\keyword{classes}
\keyword{interval data}
