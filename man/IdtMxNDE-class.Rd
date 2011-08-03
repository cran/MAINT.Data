\name{IdtMxNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMxNDE-class}

\title{Class "IdtMxNDE"}
\description{IdtMxNDE contains the results of a mixture model estimation for the Normal model, with the five different possible variance-covariance configurations.
}

\section{Slots}{
  \describe{
    \item{\code{Hmcdt}:}{Indicates whether we consider an homocedastic (TRUE) or a hetereocedasic model (FALSE)}
    \item{\code{mleNmuE}:}{Matrix with the maximum likelihood mean vectors estimates by group (each row refers to a group)}
    \item{\code{mleNmuEse}:}{Matrix with the maximum likelihood means' standard errors  by group (each row refers to a group)}
    \item{\code{Configurations}:}{List of the considered configurations}
\item{\code{grouping}:}{Factor indicating the group to which each observation belongs to}
    \item{\code{ModelNames}:}{The model acronym, indicating the model type, currently  N for Normal 
and the configuration (C1 to C5)}
\item{\code{ModelType}:}{Indicates the model;  currently, only Gaussian (Normal) distributions are implemented}
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
Class \code{"\linkS4class{IdtMxE}"}, directly.
Class \code{"\linkS4class{IdtNDE}"}, directly.
Class \code{"\linkS4class{IdtE}"}, by class "IdtMxE", distance 2.
}

\references{
Brito, P., Duarte Silva, A. P.  (2011): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics (forthcoming).
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\seealso{
  \code{\link{MANOVA}}, \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}
