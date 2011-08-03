\name{Idtqda-class}
\Rdversion{1.1}
\docType{class}
\alias{Idtqda-class}
\alias{predict,Idtqda-method}
\alias{show,Idtqda-method}

\title{Class "Idtqda"}
\description{Idtqda contains the results of Quadratic Discriminant Analysis for the interval data
}

\section{Slots}{
  \describe{
    \item{\code{prior}:}{Prior probabilities of class membership}
    \item{\code{means}:}{Matrix with the mean vectors for each group}
    \item{\code{scaling}:}{For each group i, scaling[,,i] is an array which transforms interval-valued observations so that within-groups covariance matrix is spherical.}
    \item{\code{ldet}:}{Vector of half log determinants of the dispersion matrix.}
    \item{\code{lev}:}{Levels of the grouping factor}
  }
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Idtqda")}: Classifies interval-valued observations in conjunction with qda.  }
    \item{show}{\code{signature(object = "Idtqda")}: show S4 method for the Idtqda-class
 }
	 }
}
\references{
Brito, P., Duarte Silva, A. P.  (2011): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics (forthcoming).\cr
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\seealso{
  \code{\link{MANOVA}}, \code{\link{lda}}, 
\code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}
