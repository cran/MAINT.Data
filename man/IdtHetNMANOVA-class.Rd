\name{IdtHetNMANOVA-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtHetNMANOVA-class}

\title{Class "IdtHetNMANOVA"}
\description{
 IdtHetNMANOVA contains the results of MANOVA tests on the interval-valued data, for the Heterocedastic Normal model.
}

\section{Slots}{
  \describe{
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{grouping}:}{Factor indicating the group to which each observation belongs to}
    \item{\code{H0res}:}{Model estimates under the null hypothesis}
    \item{\code{H1res}:}{Model estimates under the alternative hypothesis}
    \item{\code{QuiSq}:}{Value of the Qui-Square statistics corresponding to the performed test}
    \item{\code{df}:}{Degrees of freedom of the Qui-Square statistics}
    \item{\code{pvalue}:}{p-value of the Qui-Square statistics value, obtained from the Qui-Square distribution with df degrees of freedom}
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
    \item{qda}{\code{signature(x = "IdtHetNMANOVA")}: performs Quadratic Discriminant Analysis using the variance-covariance matrices estimated under Gaussian heterocedastic models.}
	 }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.}

\seealso{
  \code{\link{MANOVA}}, \code{\link{qda}}, 
\code{\linkS4class{IData}} 
}
\keyword{classes}
\keyword{interval data}
\keyword{MANOVA}


