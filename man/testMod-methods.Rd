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
In the former case a Restricted model is always compared against the most parcimonious alternative that encompasses it, and in latter all
possible comparisons are performed}

}

\value{An object of type ConfTests with the results of the tests performed}

\keyword{methods}
\keyword{interval data}
\keyword{likelihood ratio test}

\seealso{ \code{\linkS4class{ConfTests}} }

\examples{
# Create an Interval-Data object containing the intervals of temperatures by quarter for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

#Estimate parameters by maximum likelihood
ChinaE <- mle(ChinaT)
cat("China maximum likelhiood estimation results =\n") ; print(ChinaE)

# Perform Likelihood-Ratio tests comparing models with consecutive nested configuration 
testMod(ChinaE)

# Perform Likelihood-Ratio tests comparing all possible models 
testMod(ChinaE,FullMod="All")

# Compare model with covariance Configuration 2 against model with covariance Configuration 1 (unrestricted covariance)  
testMod(ChinaE,RestMod=2,FullMod=1)

}


