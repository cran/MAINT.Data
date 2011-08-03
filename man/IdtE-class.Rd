\name{IdtE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtE-class}
\alias{BestModel,IdtE-method}
\alias{show,IdtE-method}
\alias{testMod,IdtE-method}

\title{Class "IdtE"}
\description{IdtE contains the results of model estimation.
}

\section{Slots}{
  \describe{
    \item{\code{ModelNames}:}{The model acronym, indicating the model type, currently N for Normal 
and the configuration (C1 to C5)}
    \item{\code{ModelType}:}{Indicates the model; currently only Gaussian (Normal) distributions are implemented}
    \item{\code{ModelConfig}:}{Configuration of the variance-covariance matrix: a vector of integers in the range 1 to 5}
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{SelCrit}:}{The model selection criterion; currently,  AIC and BIC are implemented}
    \item{\code{logLiks}:}{The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Value of the AIC criterion }
    \item{\code{BICs}:}{Value of the BIC criterion }
    \item{\code{BestModel}:}{The index of best model according to the chosen selection criterion}
  }
}

\section{Methods}{
  \describe{
    \item{BestModel}{\code{signature(Idt = "IdtE")}: Selects the best model according to the chosen selection criterion (currently, AIC or BIC) }
    \item{show}{\code{signature(object = "IdtE")}: show S4 method for the IDtE-class
 }
  \item{testMod}{\code{signature(Idt = "IdtE")}: Performs statistical likelihood-ratio tests that evaluate the goodness-of-fit of a nested model against a more general one.
 }
	 }
}

\references{
Brito, P., Duarte Silva, A. P.  (2011): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics (forthcoming).
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\seealso{
  \code{\link{mle}}, \code{\link{MANOVA}}, \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}

