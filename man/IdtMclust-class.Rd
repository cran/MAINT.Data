\name{IdtMclust-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMclust-class}
\alias{show,IdtMclust-method}

\title{Class IdtMclust}
\description{IdtMclust contains the results of fitting mixtures of Gaussian distributions to interval data represented by objects of class \code{\linkS4class{IData}}.}

\section{Slots}{
  \describe{
    \item{\code{call}:}{The matched call that created the IdtMclust object} 
    \item{\code{data}:}{The IData data object} 
    \item{\code{NObs}:}{Number of entities under analysis (cases)}
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{SelCrit}:}{The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{Hmcdt}:}{Indicates whether the optimal model corresponds to a homocedastic (TRUE) or a hetereocedasic (FALSE) setup}
    \item{\code{BestG:}}{The optimal number of mixture components.}
    \item{\code{BestC:}}{The configuration case of the variance-covariance matrix in the optimal model}
    \item{\code{logLiks}:}{The logarithms of the likelihood function for the different models tried}
    \item{\code{logLik}:}{The logarithm of the likelihood function for the optimal model}
    \item{\code{AICs}:}{The values of the AIC criterion for the different models tried}
    \item{\code{aic}:}{The value of the AIC criterion for the he optimal model}
    \item{\code{BICs}:}{The values of the BIC criterion for the different models tried}
    \item{\code{bic}:}{The value of the BIC criterion for the he optimal model}
    \item{\code{parameters}}{
      A list with the following components:
      \describe{
        \item{pro}{A vector whose \emph{k}th component is the mixing proportion for the \emph{k}th component of the mixture model.}
        \item{mean}{The mean for each component. If there is more than one component, this is a matrix whose kth column is the mean of the \emph{k}th component of the mixture model.}
        \item{covariance}{A three-dimensional array with the covariance estimates. If Hmcdt is FALSE (heterocedastic setups) the third dimension levels run through the BestG mixture components, with one different covariance matrix for each level. Otherwise (homocedastic setups), there is only one covariance matrix and the size of the third dimension equals one.}
      }
    }
    \item{\code{z:}}{A matrix whose \emph{[i,k]}th entry is the probability that observation \emph{i} in the test data belongs to the \emph{k}th class.}
    \item{\code{classification:}}{The classification corresponding to \code{z}, i.e. \code{map(z)}.}
    \item{\code{allres:}}{A list with the detailed results for all models fitted.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "IdtMclust")}: show S4 method for the IdtMclust-class}
    \item{summary}{\code{signature(object = "IdtMclust")}: summary S4 method for the IdtMclust-class}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Brito, P., Duarte Silva, A. P. and Dias, J. G. (2015), Probabilistic Clustering of Interval Data. \emph{Intelligent Data Analysis} \bold{19}(2), 293--313.
}

\seealso{
  \code{\link{Idtmclust}}  
}

\keyword{classes}
\keyword{interval data}
\keyword{model based clustering}

