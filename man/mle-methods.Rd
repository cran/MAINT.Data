\name{mle-methods}
\docType{methods}
\alias{mle-methods}
\alias{mle,IData-method}
\alias{mle}   
\title{ Methods for Function mle in Package `MAINT.Data'}

\description{Performs maximum likelihood estimation for different distributions}

\section{Methods}{
\describe{

\item{\code{signature(Idt = "IData")}}{
Performs maximum likelihood estimation for interval-valued data, considering different distributions, and alternative variance-covariance matrix configurations}
}}


\usage{mle(Idt,Model="Normal",Config=1:5,SelCrit=c("AIC","BIC"))}

\arguments{
  \item{Idt}{An IData object representing NObs interval-valued entities.}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges; currently, only Gaussian (Normal) distributions are implemented}

   \item{Config}{Configuration of the variance-covariance matrix: a vector of integers between 1 and 5}

  \item{SelCrit}{The model selection criterion}
}

\value{An object of type IdtSngDE, containing the estimation results}

\keyword{methods}
\keyword{maximum likelihood estimation}
\keyword{interval data}

\seealso{ \code{\linkS4class{IdtSngDE}} }

\examples{
# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate parameters by maximum likelihood, assuming that 
# one of C3, C4 or C5 restricted configurations holds
ChinaE <- mle(ChinaT,Config=3:5)
cat("China maximum likelhiood estimation results =\n") ; print(ChinaE)
cat("Standard Errors of Estimators:\n") ; print(stdEr(ChinaE))
}
