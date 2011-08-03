\name{IdtNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtNDE-class}

\title{Class "IdtNDE"}
\description{IdtNDE is a class union used to store the estimation results for Gaussian (Normal) models of Interval Data.} 

\references{
Brito, P., Duarte Silva, A. P.  (2011): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics (forthcoming).
}

\seealso{
  \code{\linkS4class{IData}}, \code{\link{mle}}, \code{\link{MANOVA}},   \code{\linkS4class{IdtSngNDE}}, \code{\linkS4class{IdtMxNDE}}  
}

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(coef = "IdtNDE")}: extracts parameter estimates from objects of type IdtNDE.}
    \item{stdEr}{\code{signature(x = "IdtNDE")}: extracts estimator standard errors from objects of type IdtNDE.}
 }
}

\keyword{classes}
\keyword{Interval Data}



