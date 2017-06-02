\name{IdtSNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtSNDE-class}

\title{Class "IdtSNDE"}
\description{IdtSNDE is a class union used for storing the estimation results of Skew-Normal modelizations for Interval Data.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(coef = "IdtSNDE")}: extracts parameter estimates from objects of type IdtSNDE}
    \item{stdEr}{\code{signature(x = "IdtSNDE")}: extracts standard errors from objects of type IdtSNDE}
    \item{vcov}{\code{signature(x = "IdtSNDE")}: extracts an asymptotic estimate of the variance-covariance matrix of the paramenters estimators for objects of type IdtSNDE}
    \item{mean}{\code{signature(x = "IdtSNDE")}: extracts the mean vector estimate from objects of type IdtSNDE}
    \item{var}{\code{signature(x = "IdtSNDE")}: extracts the variance-covariance matrix estimate from objects of type IdtSNDE}
    \item{cor}{\code{signature(x = "IdtSNDE")}: extracts the correlation matrix estimate from objects of type IdtSNDE}
 }
}

\seealso{
  \code{\linkS4class{IData}}, \code{\link{mle}}, \code{\link{MANOVA}},   \code{\linkS4class{IdtSngSNDE}}, \code{\linkS4class{IdtMxSNDE}}  
}

\keyword{classes}
\keyword{Interval Data}



