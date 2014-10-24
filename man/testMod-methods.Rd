\name{testMod-methods}
\docType{methods}
\alias{testMod-methods}
\alias{testMod}   
\title{ Methods for Function testMod in Package `MAINT.Data'}

\description{Performs statistical likelihood-ratio tests that evaluate the goodness-of-fit of a nested model against a more general one. }

\usage{testMod(Idt,RestMod=1:length(Idt@ModelNames),FullMod="Next")}

\arguments{
  \item{Idt}{An IData object representing NObs interval-valued entities}

  \item{RestMod}{Indices of the restricted models beeing evaluated in the NULL hypothesis}

  \item{FullMod}{Either indices of the general models beeing evaluated in the alternative hypothesis or the strings "Next" (default) or "All".
In the former case a Restricted model is always compared against the most parsimonious alternative that encompasses it, and in latter all
possible comparisons are performed}

}

\value{An object of type ConfTests with the results of the tests performed}

\keyword{methods}
\keyword{interval data}
\keyword{likelihood ratio test}

\seealso{ \code{\linkS4class{ConfTests}} }

\examples{

# Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate by maximum likelihood the parameters of Gaussian models 
# with configurations C1, C3, C4 and C5 for the Winter (1st and 4th) quarter intervals
#ChinaWTE <- mle(ChinaT[,c(1,4)],Config=c(1,3,4,5))
ChinaWTE <- mle(ChinaT[,c(1,4)])
cat("China maximum likelhiood estimation results for Winter quarters:\n")
print(ChinaWTE)

# Perform Likelihood-Ratio tests comparing models with consecutive nested Configuration 
testMod(ChinaWTE)

# Perform Likelihood-Ratio tests comparing all possible models 
testMod(ChinaWTE,FullMod="All")

# Compare model with covariance Configuration 4 
# against model with covariance Configuration 1 (unrestricted covariance)  
testMod(ChinaWTE,RestMod=4,FullMod=1)

}



