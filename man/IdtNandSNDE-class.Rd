\name{IdtNandSNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtNandSNDE-class}

\title{Class "IdtNandSNDE"}
\description{IdtNandSNDE is a class union used for storing the estimation results of Normal and Skew-Normal modelizations for Interval Data.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(coef = "IdtNandSNDE")}: extracts parameter estimates from objects of type IdtNandSNDE.}
    \item{stdEr}{\code{signature(x = "IdtNandSNDE")}: extracts standard errors from objects of type IdtNandSNDE.}
    \item{vcov}{\code{signature(x = "IdtNandSNDE")}: extracts an estimate of the variance-covariance matrix of the parameters estimators for objects of type IdtNandSNDE.}
 }
}

\seealso{
  \code{\linkS4class{IData}}, \code{\link{mle}}, \code{\link{fasttle}}, \code{\link{fulltle}}, 
  \code{\link{MANOVA}}, \code{\link{RobMxtDEst}}, \code{\linkS4class{IdtSngNandSNDE}}, \code{\linkS4class{IdtMxNandSNDE}}  
}

\keyword{classes}
\keyword{Interval Data}



