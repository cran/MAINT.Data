\name{IdtNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtNDE-class}

\title{Class "IdtNDE"}
\description{IdtNDE is a class union used for storing the estimation results of Normal modelizations for Interval Data.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
  \code{\linkS4class{IData}}, \code{\link{mle}}, \code{\link{fasttle}}, \code{\link{fulltle}}, 
  \code{\link{MANOVA}}, \code{\link{RobMxtDEst}},  \code{\linkS4class{IdtSngNDE}}, \code{\linkS4class{IdtMxNDE}}
}  

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(coef = "IdtNDE")}: extracts parameter estimates from objects of type IdtNDE}
    \item{stdEr}{\code{signature(x = "IdtNDE")}: extracts standard errors from objects of type IdtNDE}
    \item{vcov}{\code{signature(x = "IdtNDE")}: extracts an estimate of the variance-covariance matrix of the paramenters estimators for objects of type IdtNDE}
 }
}

\keyword{classes}
\keyword{Interval Data}



