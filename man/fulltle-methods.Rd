\name{fulltle-methods}
\docType{methods}
\alias{fulltle-methods}
\alias{fulltle,IData-method}
\alias{fulltle}   % Check if this sintax is really correct !!!
\title{ Methods for Function fulltle in Package \sQuote{MAINT.Data}}

\description{Performs maximum trimmed likelihood estimation by an exact algorithm (full enumeratiom of all k-trimmed subsets)}

\section{Methods}{
\describe{

\item{\code{signature(Idt = "IData")}}{
Performs maximum trimmed likelihood estimation for interval-valued data using an exact algorithm, and assuming a Gaussian distribution, and considering alternative variance-covariance matrix configurations}
}}

\keyword{methods}
\keyword{trimmed maximum likelihood estimator}
\keyword{interval data}

\usage{fulltle(Idt, alpha=0.75, reweighted=TRUE, CorrF=c("smallsmp","consistent","none"),
  outlin=c("MidPandLogR","MidP","LogR"),
  CovCase=1:4, SelCrit=c("BIC","AIC"),
  force=FALSE,otpType=c("OnlyEst","SetMD2andEst"), \dots)
}

\arguments{
    \item{Idt}{An IData object representing interval-valued entities.}

    \item{alpha}{Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*Idt@NIVar observations are used for computing the trimmed likelihood. Allowed values are between 0.5 and 1.}

    \item{reweighted}{Should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}

    \item{CorrF}{Whether to use sample correction factors for the covariance estimates; options are \dQuote{smallsmp} (default), for small sample correction, \dQuote{consistent}, for consistency correction and \dQuote{none} for no correction.}

    \item{outlin}{The type of outliers to be consideres. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}

    \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

    \item{SelCrit}{The model selection criterion.}

    \item{force}{A boolean flag indicating whether, for moderate or large data sets the algorithm should proceed anyway, regardless of an expected long excution time, due to exponential explosions in the number of different subsets that need to be avaluated by fulltle}

    \item{otpType}{The amount of output returned by fulltle. Current options are \dQuote{OnlyEst} (default) where only an \sQuote{IdtE} object with the fulltle estimates is returned and \dQuote{SetMD2andEst} which returns a list with an \sQuote{IdtE} object of fulltle estimates, a vector with the final trimmed subset elements used to compute these estimates and the corresponding robust squared Mahalanobis distances.}

    \item{\dots}{Further arguments to be passed to internal functions of \sQuote{fulltle}.} 
}

\value{

If argument \sQuote{otpType} is set to \dQuote{OnlyEst}, an object of class \sQuote{IdtE} with the fulltle estimates, their log-likelihood values, and the value of the comparison criterion used to select the covariance configurations.

If argument \sQuote{otpType} is set to \dQuote{SetMD2andEst} a list with the following components:

\item{sol}{An object of class \sQuote{IdtE} with the fulltle estimates, their log-likelihood values, and the value of the comparison criterion used to select the covariance configurations.} 
\item{Set}{A vector with the final trimmed subset elements used to compute the fulltle estimates.} 
\item{RobMD2}{A vector with the robust squared Mahalanobis distances used to select the trimmed subset.} 

} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Hadi, A. S. and Luceno, A. (1997), Maximum trimmed likelihood estimators: a unified approach, examples, and algorithms.
\emph{Computational Statistics and Data Analysis} \bold{25}(3), 251--272.\cr
}


\seealso{ \code{\link{fasttle}} }


\examples{

# Create an Interval-Data object containing the intervals for characteristics 
# of 27 cars models.

Cars <- IData(CarsData[1:8],Seq="MidPLogR_VarbyVar",
  VarNames=c("Price","EngineCapacity","TopSpeed","Acceleration"))

#Display the first and last observations

# Estimate normal distributuion  parameters robustly by full maximum trimmed likelihood, 

CarsTE <- fulltle(Cars)
cat("Cars data -- normal maximum trimmed likelhiood estimation results:\n")
print(CarsTE)

}
