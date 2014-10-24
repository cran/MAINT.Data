\name{IdtMxE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMxE-class}

\title{Class "IdtMxE"}
\description{IdtMxE contains the results of a mixture model estimation. 
}

\section{Slots}{
  \describe{
    \item{\code{grouping}:}{Factor indicating the group to which each observation belongs to}
    \item{\code{ModelNames}:}{The model acronym, indicating the model type, currently N for Normal 
followed by the configuration (C1 to C5)}
    \item{\code{ModelType}:}{Indicates the model; currently, only Gaussian (Normal) distributions are implemented}
    \item{\code{ModelConfig}:}{Configuration of the variance-covariance matrix: a vector of integers in the range 1 to 5}
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{SelCrit}:}{The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{logLiks}:}{The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Value of the AIC criterion }
    \item{\code{BICs}:}{Value of the BIC criterion }
    \item{\code{BestModel}:}{The index of best model according to the chosen selection criterion}
  }
}
\section{Extends}{
Class \code{"\linkS4class{IdtE}"}, directly.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.}

\seealso{
  \code{\link{MANOVA}}, \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}

