\name{IdtMxNDRE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMxNDRE-class}

\title{Class IdtMxNDE}
\description{IdtMxNDRE contains the results of a mixture Normal model robust parameter estimation, with the four different possible variance-covariance configurations.}

\section{Slots}{
  \describe{
    \item{\code{Hmcdt}:}{Indicates whether we consider an homocedastic (TRUE) or a hetereocedasic model (FALSE)}
    \item{\code{RobNmuE}:}{Matrix with the robust mean vectors estimates by group (each row refers to a group)}
    \item{\code{CovConfCases}:}{List of the considered configurations}
    \item{\code{grouping}:}{Inherited from class \code{\linkS4class{IdtMxE}}. Factor indicating the group to which each observation belongs to}
    \item{\code{ModelNames}:}{Inherited from class \code{\linkS4class{IdtE}}. The model acronym formed by a "N", indicating a Normal model, followed by the configuration (Case 1 through Case 4)}
    \item{\code{ModelType}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the model; always set to "Normal" in objects of the IdtMxNDRE class}    
\item{\code{ModelConfig}:}{Inherited from class \code{\linkS4class{IdtE}}. Configuration case of the variance-covariance matrix: Case 1 through Case 4}    
    \item{\code{NIVar}:}{Inherited from class \code{\linkS4class{IdtE}}. Number of interval variables}
    \item{\code{SelCrit}:}{Inherited from class \code{\linkS4class{IdtE}}. The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{logLiks}:}{Inherited from class \code{\linkS4class{IdtE}}. The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the AIC criterion }
    \item{\code{BICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the BIC criterion }
    \item{\code{BestModel}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the best model according to the chosen selection criterion}
    \item{\code{SngD}:}{Inherited from class \code{\linkS4class{IdtE}}. Boolean flag indicating whether a single or a mixture of distribution were estimated. Always set to FALSE in objects of class IdtMxNDRE }
    \item{\code{Ngrps}:}{Inherited from class \code{\linkS4class{IdtMxE}}. Number of mixture components}
  }
}  

\section{Extends}{
Class \code{\linkS4class{IdtMxE}}, directly.
Class \code{\linkS4class{IdtE}}, by class \code{\linkS4class{IdtMxE}}, distance 2.
}

\section{Methods}{
No methods defined with class IdtMxNDRE in the signature.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Hadi, A. S. and Luceno, A. (1997), Maximum trimmed likelihood estimators: a unified approach, examples, and algorithms.
\emph{Computational Statistics and Data Analysis} \bold{25}(3), 251--272.
}

\seealso{
  \code{\linkS4class{IdtE}}, \code{\linkS4class{IdtMxE}}, \code{\linkS4class{IdtMxNDE}}, \code{\linkS4class{IdtSngNDRE}}, \code{\link{RobMxtDEst}}, \code{\linkS4class{IData}}
}

\keyword{classes}
\keyword{interval data}