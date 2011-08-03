\name{MANOVA-methods}
\docType{methods}
\alias{MANOVA}   % Check if this sintax is really correct !!!
\alias{MANOVA-methods}
\alias{MANOVA,IData-method}
\title{Methods for Function MANOVA in Package `MAINT.Data'}
\description{Performs MANOVA tests based on likelihood ratios allowing for both homocedastic and heterocedastic setups
}
\section{Methods}{
\describe{

\item{\code{signature(Idt = "IData")}}{
Performs MANOVA tests on interval-valued data}
}}

\keyword{methods}
\keyword{MANOVA}
\keyword{interval data}
\keyword{likelihood ratio tets}

\usage{MANOVA(Idt, grouping, Model="Normal", Config=1:5,
			SelCrit=c("AIC","BIC"), Mxt=c("Hom","Het"), tol=1.0e-4)}

\arguments{
  \item{Idt}{An IData object representing NObs interval-valued entities}

  \item{grouping}{Factor indicating the group to which each observation belongs to}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges; currently, only Gaussian (Normal) distributions are implemented}

   \item{Config}{Configuration of the variance-covariance matrix: a vector of integers between 1 and 5}

   \item{SelCrit}{The model selection criterion}

   \item{Mxt}{Indicates whether we consider an homocedastic ("Hom") or a hetereocedasic ("Het") model configuration  }

   \item{tol}{Tolerance level for the standard deviation of non-constant variables. When a MidPoint or LogRange has a within-groups standard deviation below tol, it is considered to be a constant}
}

\value{An object of type IdtMANOVA, containing the estimation and test results}

\seealso{ \code{\linkS4class{IdtMANOVA}} }

\examples{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.
ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#MANOVA tests, assuming that one of C3, C4 or C5 restricted configurations holds

#Classical (homocedastic) MANOVA tests
ManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg,Config=3:5)
cat("China, MANOVA by geografical regions results =\n") 
print(ManvChina)

#Heterocedastic MANOVA tests
HetManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg,Mxt="Het",Config=3:5)
cat("China, heterocedastic MANOVA by geografical regions results =\n") 
print(HetManvChina)
}
