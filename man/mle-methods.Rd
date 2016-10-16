\name{mle-methods}
\docType{methods}
\alias{mle-methods}
\alias{mle,IData-method}
\alias{mle}   
\title{ Methods for Function mle in Package `MAINT.Data'}

\description{Performs maximum likelihood estimation for parametric models of interval data}

\usage{

\S4method{mle}{IData}(Idt, Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4,
  SelCrit=c("BIC","AIC"), OptCntrl=list(), \dots)

}

\arguments{
  \item{Idt}{An IData object representing interval-valued entities.}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges. Current alternatives are \dQuote{Normal} for Gaussian,
distributions, \dQuote{SNNormal} for Skew-Normal and \dQuote{NrmandSKN} for both Gaussian and Skew-Normal distributions.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

   \item{SelCrit}{The model selection criterion.}

   \item{OptCntrl}{List of optional control parameters to passed to the optimization routine. See the documentation of RepLOptim for a description of the available options.}

  \item{\dots}{Other named arguments.}
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.\cr
}

\seealso{
\code{\link{IData}}, \code{\link{RepLOptim}}
}

\keyword{methods}
\keyword{maximum likelihood estimation}
\keyword{interval data}

\examples{
# Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate parameters by maximum likelihood, assuming that one of 
# the C2, C3 or C4 restricted Covariance configurations holds

ChinaE <- mle(ChinaT,CovCase=2:4)
cat("China maximum likelhiood estimation results =\n")
print(ChinaE)
cat("Standard Errors of Estimators:\n")
print(stdEr(ChinaE))
}
