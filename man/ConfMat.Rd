\name{ConfMat}

\title{Confussion Matrices for classification results}

\alias{ConfMat}

\description{
\sQuote{ConfMat} creates confussion matrices from two factor describing, respectively, original classes and predicted  classification results }

\usage{ConfMat(origcl, predcl, otp=c("absandrel","abs","rel"), dec=3)}


\arguments{

  \item{origcl}{A factor describing the original classes.}

  \item{predcl}{A factor describing the predicted classes.}

  \item{otp}{A string describing the output to be displayed and returned. Alternatives are \dQuote{absandrel} for two confusion matrices, respectively with absolute and relative frequencies, \dQuote{abs} for a confusion matrix with absolute frequencies, and \dQuote{rel} for a confusion matrix relative frequencies.}  

  \item{dec}{The number of decimal digits to display in matrices of relative frequencies.} 
}

\value{When argument \sQuote{otp} is set to \dQuote{absandrel} (default), a list with two confusion matrices, respectively with absolute and relative frequencies. When argument \sQuote{otp} is set to \dQuote{abs} a confusion matrix with absolute frequencies, and when argument \sQuote{otp} is set to \dQuote{rel} a confusion matrix with relative frequencies.}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{lda}}, \code{\link{qda}}, \code{\link{snda}}, \code{\link{Roblda}}, \code{\link{Robqda}}, \code{\link{DACrossVal}}}

\examples{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#Linear Discriminant Analysis

ChinaT.lda <- lda(ChinaT,ChinaTemp$GeoReg)
ldapred <- predict(ChinaT.lda,ChinaT)$class

# lda resubstitution confusion matrix

ConfMat(ChinaTemp$GeoReg,ldapred)

#Quadratic Discriminant Analysis

ChinaT.qda <- qda(ChinaT,ChinaTemp$GeoReg)
qdapred <- predict(ChinaT.qda,ChinaT)$class

# qda resubstitution confusion matrix

ConfMat(ChinaTemp$GeoReg,qdapred)


}

\keyword{discriminant analysis} 
\keyword{error rates} 

