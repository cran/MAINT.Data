\name{AgrMcDt}
\alias{AgrMcDt}
\title{Agregate Micro Data}
\description{
AgrMcDt creates \code{\linkS4class{IData}} objects by agregating a Data Frame of Micro Data.
}

\usage{
AgrMcDt(MicDtDF, agrby, agrcrt="minmax")
}

\arguments{
  \item{MicDtDF}{A data frame with the original values of the micro data.}
  \item{agrby}{A factor with categories on which the micro data should be agregated.}
  \item{agrcrt}{The agregation criterion. Either the \sQuote{minmax} string, or a two dimensional vecto with the prob. value for the left (lower) percentile, followed by the prob. value for the right (upper) percentile, used in the agregation.}
}

\value{An object of class \code{\linkS4class{IData}} with the data set of Interval-valued variables resulting from the agregation performed.}

\seealso{
  \code{\linkS4class{IData}} 
}

\keyword{interval data}
