\name{IdtClMANOVA-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtClMANOVA-class}

\title{Class "IdtClMANOVA"}
\description{
  IdtClMANOVA contains the results of classical MANOVA tests applied to the chosen variables
}

\section{Slots}{
  \describe{
    \item{\code{NIVar}:}{Number of interval variables }
  \item{\code{grouping}:}{Factor indicating the group to which the observation belongs to}
    \item{\code{H0res}:}{Model estimates under the null hypothesis}
    \item{\code{H1res}:}{Model estimates under the alternative hypothesis}
    \item{\code{QuiSq}:}{Value of the Qui-Square statistics corresponding to the performed test}
    \item{\code{df}:}{Degrees of freedom of the Qui-Square statistics}
    \item{\code{pvalue}:}{p-value of the Qui-Square statistics value, 
obtained from the Qui-Square distribution with df degrees of freedom}
    \item{\code{H0logLik}:}{Logarithm of the Likelihood function under the null hypothesis}
    \item{\code{H1logLik}:}{Logarithm of the Likelihood function under the alternative hypothesis}
  }
}

\section{Extends}{
Class \code{"\linkS4class{IdtMANOVA}"}, directly.
Class \code{"\linkS4class{LRTest}"}, by class "IdtMANOVA", distance 2.
}

\section{Methods}{
  \describe{
    \item{lda}{\code{signature(x = "IdtClMANOVA")}: lda performs Linear Discriminat Analysis using the estimated model parameters}
	 }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.}

\seealso{
  \code{\link{MANOVA}}, \code{\link{lda}}, 
\code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}
\keyword{MANOVA}


