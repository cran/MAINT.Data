\name{IData}
\alias{IData}
\alias{is.IData}
\title{Interval Data objects}
\description{
\sQuote{IData} creates IData objects from data frames of interval bounds or MidPoint/LogRange values of the interval-valued observations.

\sQuote{is.IData} tests if its argument is an IData object.
}

\usage{
IData(Data, 
Seq = c("LbUb_VarbyVar", "MidPLogR_VarbyVar", "AllLb_AllUb",
"AllMidP_AllLogR"), VarNames=NULL, ObsNames=row.names(Data))

is.IData(object)
}

\arguments{
  \item{Data}{a data frame of interval bounds or MidPoint/LogRange values.}

  \item{Seq}{the format of \sQuote{Data} data frame. Available options are:\cr
\dQuote{LbUb_VarbyVar}: lower bounds followed by upper bounds, variable by variable.\cr
\dQuote{MidPLogR_VarbyVar}: MidPoints followed by LogRanges, variable by variable.\cr
\dQuote{AllLb_AllUb}: all lower bounds followed by all upper bounds, in the same variable order.\cr
\dQuote{AllMidP_AllLogR}: all MidPoints followed all LogRanges, in the same variable order.\cr
}

  \item{VarNames}{An optional vector of names to be assigned to the Interval-Valued Variables.}

  \item{ObsNames}{An optional vector of names assigned to the individual observations.}

  \item{object}{An R object}
}

\details{Objects of type \sQuote{IData}, describe a data set of \sQuote{NObs} observations on \sQuote{NIVar} Interval-Valued variables. This function creates an interval-data object from a data-frame with either the lower and upper bounds of the observed intervals or by their midpoints and log-ranges.}

\keyword{interval data}

\examples{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

cat("Summary of the ChinaT IData object:\n")  ; print(summary(ChinaT))
cat("ChinaT first ant last observations:\n")  
print(head(ChinaT,n=3))
cat("\n...\n")
print(tail(ChinaT,n=3))
}



