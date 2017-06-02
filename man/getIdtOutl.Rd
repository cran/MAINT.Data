\name{getIdtOutl}
\alias{getIdtOutl}

\title{Get Interval Data Outliers}

\description{Identifies outliers in a data set of Interval-valued variables}

\usage{getIdtOutl(Idt, IdtE=NULL, muE=NULL, SigE=NULL,
  eta=0.025, Rewind=IdtE$RewghtdSet, m=length(Rewind),
  RefDist=c("ChiSq","HardRockeAdjF","HardRockeAsF","CerioliBetaF"),
  multiCmpCor=c("never","always","iterstep"), 
  outlin=c("MidPandLogR","MidP","LogR"))
}


\arguments{
    \item{Idt}{An IData object representing interval-valued entities.}
    \item{IdtE}{Ao object of class \sQuote{IdtSngNDRE} or \sQuote{IdtSngNDE} containing mean and covariance estimates.}
    \item{muE}{Vector with the mean estimates used to find Mahalanobis distances.}
    \item{SigE}{Matrix with the covariance estimates used to find Mahalanobis distances.}
    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier.}
    \item{Rewind}{A vector with the subset of entities used to compute trimmed mean and covariance estimates. Only used when the \sQuote{RefDist} argument is set to \dQuote{CerioliBetaF.}}
    \item{m}{Number of entities used to compute trimmed mean and covariance estimates. Not used when the \sQuote{RefDist} argument is set to \dQuote{ChiSq.}}
    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at \sQuote{eta}. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n); and \sQuote{iterstep} -- as sugested by Cerioli (2010), make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n), and if no outliers were detected stop. Otherwise, make a second step testing for outliers at \sQuote{eta}. }
    \item{RefDist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq},\dQuote{HardRockeAsF}, \dQuote{HardRockeAdjF} and \dQuote{CerioliBetaF}, respectivelly for the usual Qui-squared, the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005), and the Beta and F distributions proposed by Cerioli (2010).}    
    \item{outlin}{The type of outliers to be considered. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}
}

\value{A vector with the indices of the entities identified as outliers.} 

\references{
Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156. 

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr
}

\seealso{ \code{\link{fasttle}}, \code{\link{fulltle} }


}
