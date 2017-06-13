\name{IdtSngNDRE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtSngNDRE-class}

\title{Class IdtSngNDRE}
\description{Contains the results of a single class robust estimation for the Normal distribution, with the four different possible variance-covariance configurations.}

\section{Slots}{
  \describe{
    \item{\code{RobNmuE}:}{Matrix with the maximum likelihood mean vectors estimates}
    \item{\code{CovConfCases}:}{List of the considered configurations}
    \item{\code{ModelNames}:}{Inherited from class \code{\linkS4class{IdtE}}. The model acronym formed by a "N", indicating a Normal model, followed by the configuration (Case 1 through Case 4)}
    \item{\code{ModelType}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the model; always set to "Normal" in objects of the IdtSngNDRE class}    
    \item{\code{ModelConfig}:}{Inherited from class \code{\linkS4class{IdtE}}. Configuration of the variance-covariance matrix: Case 1 through Case 4}    
    \item{\code{NIVar}:}{Inherited from class \code{\linkS4class{IdtE}}. Number of interval variables}
    \item{\code{SelCrit}:}{Inherited from class \code{\linkS4class{IdtE}}. The model selection criterion; currently, AIC and BIC are implemented}    
    \item{\code{logLiks}:}{Inherited from class \code{\linkS4class{IdtE}}. The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the AIC criterion }
    \item{\code{BICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the BIC criterion }
    \item{\code{BestModel}:}{Inherited from class \code{\linkS4class{IdtE}}. Bestmodel indicates the best model according to the chosen selection criterion}
    \item{\code{SngD}:}{Inherited from class \code{\linkS4class{IdtE}}. Boolean flag indicating whether a single or a mixture of distribution were estimated. Always set to TRUE in objects of class IdtSngNDRE }
 }
}

\section{Extends}{
Class \code{\linkS4class{IdtSngDE}}, directly.
Class \code{\linkS4class{IdtE}}, by class \code{\linkS4class{IdtSngDE}}, distance 2.
}

\section{Methods}{
No methods defined with class IdtSngNDRE in the signature.
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
  \code{\linkS4class{IData}}, \code{\link{fasttle}}, \code{\link{fulltle}}, \code{\linkS4class{IdtSngNDE}},  \code{\linkS4class{IdtMxNDRE}}  
}

\keyword{classes}
\keyword{interval data}
