\name{fasttle-methods}
\docType{methods}
\alias{fasttle-methods}
\alias{fasttle,IData-method}
\alias{fasttle}
\title{ Methods for Function fasttle in Package \sQuote{MAINT.Data}}

\description{Performs maximum trimmed likelihood estimation by the fasttle algorithm}

\section{Methods}{
\describe{

\item{\code{signature(Idt = "IData")}}{
Performs maximum trimmed likelihood estimation for interval-valued data using the fasttle algorithm, assuming a Gaussian distribution, and considering alternative variance-covariance matrix configurations}
}}

\keyword{methods}
\keyword{trimmed maximum likelihood estimator}
\keyword{interval data}

\usage{fasttle(Idt,
    CovCase=1:4,
    SelCrit=c("BIC","AIC"),
    alpha=control@alpha,
    nsamp = control@nsamp,
    seed=control@seed,
    trace=control@trace,
    use.correction=control@use.correction,
    ncsteps=control@ncsteps,
    getalpha=control@getalpha,
    rawMD2Dist=control@rawMD2Dist,				
    MD2Dist=control@MD2Dist,
    eta=control@eta,
    multiCmpCor=control@multiCmpCor,				
    getkdblstar=control@getkdblstar,
    outlin=control@outlin,
    trialmethod=control@trialmethod,
    m=control@m,
    reweighted = control@reweighted,
    otpType=control@otpType,
    control=RobEstControl(), \dots)
}

\arguments{
    \item{Idt}{An IData object representing interval-valued entities.}
    \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}
    \item{SelCrit}{The model selection criterion.}
    \item{alpha}{Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*Idt@NIVar observations are used for computing the trimmed likelihood. Allowed values are between 0.5 and 1.}
    \item{nsamp}{Number of subsets used for initial estimates. Note that when argument \sQuote{getalpha} is set to \dQuote{TwoStep}
the final value of \sQuote{alpha} is estimated by a two-step procedure and the value of argument \sQuote{alpha} is only used
to specify the size of the samples used in the first step.}  
    \item{seed}{Initial seed for random generator, like \code{\link{.Random.seed}}, see \code{\link{rrcov.control}}.}
    \item{trace}{Logical (or integer) indicating if intermediate results should be printed; defaults to \code{FALSE}.}
    \item{use.correction}{ whether to use finite sample correction factors; defaults to \code{TRUE}.}
    \item{ncsteps}{The maximum number of concentration steps used each iteration of the fasttle algorithm.}
    \item{getalpha}{Argument specifying if the \sQuote{alpha} parameter (roughly the percentage of the sample used for computing the trimmed likelihood) should be estimated from the data, or if the value of the argument \sQuote{alpha} should be used instead. When set to \dQuote{TwoStep}, \sQuote{alpha} is estimated by a two-step procedure with the value of argument \sQuote{alpha} specifying the size of the samples used in the first step. Otherwise, with the value of argument \sQuote{alpha} is used directly.}   
    \item{rawMD2Dist}{The assumed reference distribution of the raw MCD squared distances, which is used to find to cutoffs defining the observations kept in one-step reweighted MCD estimates. Alternatives are \sQuote{ChiSq},\sQuote{HardRockeAsF} and \sQuote{HardRockeAdjF}, respectivelly for the usual Qui-squared, and the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005).}
    \item{MD2Dist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Qui-squared, the Beta and F distributions proposed by Cerioli (2010).}    
    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier. Defines the raw MCD Mahalanobis distances cutoff used to choose the observations kept in the reweightening step.}
    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at \sQuote{eta}. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n); and \sQuote{iterstep} -- as sugested by Cerioli (2010), make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n), and if no outliers were detected stop. Otherwise, make a second step testing for outliers at \sQuote{eta}. }   
    \item{getkdblstar}{Argument specifying the size of the initial small (in order to minimize the probability of outliers) subsets. If set to the string \dQuote{Twopplusone} (default) the initial sets have twice the number of interval-value variables plus one (i.e., they are the smaller samples that lead to a non-singular covariance estimate). Otherwise, an integer with the size of the initial sets.}
    \item{outlin}{The type of outliers to be considered. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}
    \item{trialmethod}{The method to find a trial subset used to initialize each replication of the fasttle algorithm. The current options are \dQuote{simple} (default) that simply selects \sQuote{kdblstar} observations at random, and \dQuote{Poolm} that divides the original sample into \sQuote{m} non-overlaping subsets, applies the \sQuote{simple trial} and the refinement methods to each one of them, and merges the results into a trial subset.}
    \item{m}{Number of non-overlaping subsets used by the trial method when the argument of \sQuote{trialmethod} is set to 'Poolm'.}
    \item{reweighted}{Should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}
    \item{otpType}{The amount of output returned by fasttle. Current options are \dQuote{OnlyEst} (default) where only an \sQuote{IdtE} object with the fasttle estimates is returned, \cr
\dQuote{SetMD2andEst} which returns a list with an \sQuote{IdtE} object of fasttle estimates, a vector with the final trimmed subset elements used to compute these estimates and the corresponding robust squared Mahalanobis distances, and \cr
\dQuote{SetMD2EstandPrfSt} wich returns a list with the previous three components plust a list of some performance statistics concerning the algorithm execution.}
  \item{control}{a list with estimation options - this includes those above provided in the function specification. See
\code{\link{RobEstControl}} for the defaults.  If \code{control} is supplied, the parameters from it will be used.
If parameters are passed also in the invocation statement, they will override the corresponding elements of the control object.}
  \item{\dots}{Further arguments to be passed to internal functions of \code{fasttle}.} 
}

\value{

If argument \sQuote{otpType} is set to \dQuote{OnlyEst}, an object of class \sQuote{IdtE} with the fasttle estimates, their log-likelihood values, and the value of the comparison criterion used to select the covariance configurations.

If argument \sQuote{otpType} is set to \dQuote{SetMD2andEst} a list with the following components:

\item{sol}{An object of class \sQuote{IdtE} with the fasttle estimates, their log-likelihood values, and the value of the comparison criterion used to select the covariance configurations.} 
\item{rawSet}{A vector with the trimmed subset elements used to compute the raw (not reweighted) MCD covariance estimate for the chosen configuration.} 
\item{RewghtdSet}{A vector with the final trimmed subset elements used to compute the fasttle estimates.} 
\item{RobMD2}{A vector with the robust squared Mahalanobis distances used to select the trimmed subset.} 
\item{cnp2}{A vector of length two containing the consistency correction factor and the finite sample correction factor of the final estimate of the covariance matrix.}
\item{raw.cov}{A matrix with the raw MCD estimator used to compute the robust squared Mahalanobis distances of RobMD2.} 
\item{raw.cnp2}{A vector of length two containing the consistency correction factor and the finite sample correction factor of the raw estimate of the covariance matrix.}

If argument \sQuote{otpType} is set to \dQuote{SetMD2EstandPrfSt} a list with the following components:

\item{sol}{An object of class \sQuote{IdtE} with the fasttle estimates, their log-likelihood values, and the value of the comparison criterion used to select the covariance configurations.} 
\item{rawSet}{A vector with the trimmed subset elements used to compute the raw (not reweighted) MCD covariance estimate for the chosen configuration.} 
\item{RewghtdSet}{A vector with the final trimmed subset elements used to compute the fasttle estimates.} 
\item{RobMD2}{A vector with the robust squared Mahalanobis distances used to select the trimmed subset.} 
\item{cnp2}{A vector of length two containing the consistency correction factor and the finite sample correction factor of the final estimate of the covariance matrix.}
\item{raw.cov}{A matrix with the raw MCD estimator used to compute the robust squared Mahalanobis distances of RobMD2.} 
\item{raw.cnp2}{A vector of length two containing the consistency correction factor and the finite sample correction factor of the raw estimate of the covariance matrix.}
\item{PerfSt}{A a list with the following components: \cr
    RepSteps: A list with one component by Covariance Configuration, containing a vector with the number of refinement steps performed by the fasttle algorithm by replication. \cr
    RepLogLik: A list with one component by Covariance Configuration, containing a vector with the best log-likelihood found be fasttle algorithm by replication. \cr
    StpLogLik: A list with one component by Covariance Configuration, containing a matrix with the evolution of the log-likelihoods found be fasttle algorithm by replication and refinement step. \cr
    }
} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156. 

Hadi, A. S. and Luceno, A. (1997), Maximum trimmed likelihood estimators: a unified approach, examples, and algorithms.
\emph{Computational Statistics and Data Analysis} \bold{25}(3), 251--272.\cr

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}

\seealso{ \code{\link{fulltle}} }


\examples{

\dontrun{
# Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate parameters using the fast trimmed maximum likelihood estimator, and assuming that one of 
# the C2, C3 or C4 restricted Covariance Cases holds

Chinafasttle <- fasttle(ChinaT,CovCase=2:4)
cat("China maximum trimmed likelihood estimation results =\n")
print(Chinafasttle)
}

}
