\name{Idtlda-class}
\Rdversion{1.1}
\docType{class}
\alias{Idtlda-class}
\alias{predict,Idtlda-method}
\alias{show,Idtlda-method}

\title{Class "Idtlda"}
\description{Idtlda contains the results of Linear Discriminant Analysis for the interval data
}

\section{Slots}{
  \describe{
    \item{\code{prior}:}{Prior probabilities of class membership}
    \item{\code{means}:}{Matrix with the mean vectors for each group}
    \item{\code{scaling}:}{Matrix which transforms observations to discriminant functions, normalized so that within groups covariance matrix is spherical.}
    \item{\code{N}:}{Number of observations}
  }
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Idtlda")}: Classifies interval-valued observations in conjunction with lda. }
    \item{show}{\code{signature(object = "Idtlda")}: show S4 method for the IDdtlda-class
 }
	 }
}


\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.}

\seealso{
  \code{\link{qda}},  \code{\link{MANOVA}},  
\code{\linkS4class{IData}}
}

\keyword{classes}
\keyword{interval data}

