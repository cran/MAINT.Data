\name{IdtMxNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMxNDE-class}
\alias{IdtMxtNDE-class}

\title{Class "IdtMxNDE"}
\description{"IdtMxNDE" contains the results of a mixture Normal model maximum likelihood parameter estimation, with the four different possible variance-covariance configurations. \cr
"IdtMxtNDE" is an union of classes "IdtMxNDE" and "IdtMxNDRE", the later one containing the results of mixture Normal model parameter estimation by robust methods. }

\section{Slots}{
  \describe{
    \item{\code{Hmcdt}:}{Indicates whether we consider an homocedastic (TRUE) or a hetereocedasic model (FALSE)}
    \item{\code{mleNmuE}:}{Matrix with the maximum likelihood mean vectors estimates by group (each row refers to a group)}
    \item{\code{mleNmuEse}:}{Matrix with the maximum likelihood means' standard errors  by group (each row refers to a group)}
    \item{\code{CovConfCases}:}{List of the considered configurations}
    \item{\code{grouping}:}{Inherited from class \code{"IdtMxE"}. Factor indicating the group to which each observation belongs to}
    \item{\code{ModelNames}:}{Inherited from class \code{"IdtE"}. The model acronym formed by a "N", indicating a Normal model, followed by the configuration (Case 1 through Case 4)}
    \item{\code{ModelType}:}{Inherited from class \code{"IdtE"}. Indicates the model; always set to "Normal" in objects of the IdtMxNDE class}    
\item{\code{ModelConfig}:}{Inherited from class \code{"IdtE"}. Configuration case of the variance-covariance matrix: Case 1 through Case 4}    
    \item{\code{NIVar}:}{Inherited from class \code{"IdtE"}. Number of interval variables}
    \item{\code{SelCrit}:}{Inherited from class \code{"IdtE"}. The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{logLiks}:}{Inherited from class \code{"IdtE"}. The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Inherited from class \code{"IdtE"}. Value of the AIC criterion }
    \item{\code{BICs}:}{Inherited from class \code{"IdtE"}. Value of the BIC criterion }
    \item{\code{BestModel}:}{Inherited from class \code{"IdtE"}. Indicates the best model according to the chosen selection criterion}
    \item{\code{SngD}:}{Inherited from class \code{"IdtE"}. Boolean flag indicating whether a single or a mixture of distribution were estimated. Always set to FALSE in objects of class \code{"IdtMxNDE"} }
    \item{\code{Ngrps}:}{Inherited from class \code{"IdtMxE"}. Number of mixture components}
  }
}  

\section{Extends}{
Class \code{"\linkS4class{IdtMxE}"}, directly.
Class \code{"\linkS4class{IdtNDE}"}, directly.
Class \code{"\linkS4class{IdtE}"}, by class "IdtMxE", distance 2.
}

\section{Methods}{
  \describe{
    \item{lda}{\code{signature(x = "IdtMxtNDE")}: Linear Discriminant Analysis using the estimated model parameters.}
    \item{qda}{\code{signature(x = "IdtMxtNDE")}: Quadratic Discriminant Analysis using the estimated model parameters.}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
  \code{\link{MANOVA}}, \code{\linkS4class{IData}}, \code{\linkS4class{IdtMxNDRE}} 
}

\keyword{classes}
\keyword{interval data}
