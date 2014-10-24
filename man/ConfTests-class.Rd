\name{ConfTests-class}
\Rdversion{1.1}
\docType{class}
\alias{ConfTests-class}
\alias{show,ConfTests-method}

\title{Class "Configuration Tests"}
\description{ConfTests contains a list of the results of statistical likelihood-ratio tests 
that evaluate the goodness-of-fit of variance-covariance matrix configurations against more general ones.
}

\section{Slots}{
  \describe{
    \item{\code{TestRes}:}{List of test results; each element is an object of type LRTest, with the following components: 
value of the Qui-Square statistics, degrees of freedom of the Qui-Square statistics, 
p-value of the Qui-Square statistics value,
logarithm of the Likelihood function under the null hypothesis, 
logarithm of the Likelihood function under the alternative hypothesis}
    \item{\code{RestModels}:}{An acronym describing the restrited model (corresponding to the null hypothesis)}
    \item{\code{FullModels}:}{An acronym describing the full model (corresponding to the alternative hypothesis)}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ConfTests")}: show S4 method for the ConfTests-class}
	 }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.\cr

Noirhomme-Fraiture, M., Brito, P. (2011): "Far Beyond the Classical Data Models: Symbolic Data Analysis". Statistical Analysis and Data Mining, Volume 4, Issue 2, 157-170.
}

\seealso{
  \code{\link{mle}}, \code{\linkS4class{IData}} 
\code{\linkS4class{LRTest}}
}

\keyword{classes}
\keyword{interval data}
