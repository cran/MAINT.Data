\name{lda-methods}
\docType{methods}
\alias{lda}
\alias{lda-methods}
\alias{lda,IData-method}
\alias{lda,IdtClMANOVA-method}
\title{Methods for Function lda in Package `MAINT.Data'}
\description{Linear Discriminant Analysis}

\section{Methods}{
\describe{
\item{\code{signature(Idt = "IData")}}{Performs Linear Discriminant Analysis on interval-valued data}
}
}

\keyword{methods}
\keyword{Linear Discrimiant Analysis}
\keyword{interval data}

\usage{

\S4method{lda}{IdtClMANOVA}(x,prior="proportions",selmodel=BestModel(H1res(x)),egvtol=1.0e-9,\dots)

\S4method{lda}{IData}(x, grouping, prior="proportions", tol=1.0e-4, 
subset=1:nrow(x) ,Config=1:5,SelCrit=c("AIC","BIC"))
}

\arguments{

  \item{x}{An IdtClMANOVA object representing the results of a (homocedastic normal-based) MANOVA of Interval Data, or an IData object representing NObs interval-valued entities}

  \item{grouping}{Factor indicating the group to which each observation belongs to}

  \item{prior}{The prior probabilities of class membership.  If unspecified,
          the class proportions for the training set are used.  If
          present, the probabilities should be specified in the order
          of the factor levels}

  \item{selmodel}{The index of model of IdtClMANOVA object to be used in the construction of the discriminant functions}

  \item{tol}{A tolerance to decide if a matrix is singular; it will reject
          variables and linear combinations of unit-variance variables
          whose standard deviation is less than tol}

  \item{egvtol}{A tolerance for positive eigenvalues of the scatter matris; 
	  eigenvalues below egvtol will be considered equal to 0, and the corresponding
	  linear discriminant functions will be ignored}

  \item{subset}{An index vector specifying the cases to be used in the
          training sample.  (NOTE: If given, this argument must be
          named}

   \item{Config}{Configuration of the variance-covariance matrix: a vector of integers between 1 and 5}

   \item{SelCrit}{The model selection criterion}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{An object of type Idtlda, containing the linear discriminant functions}

\seealso{ \code{\linkS4class{Idtlda}}, \code{\link{DACrossVal}}, \code{\link{qda}} }

\examples{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.
ChinaT <- IData(ChinaTemp[1:8])

#Linear Discriminant Analysis, assuming that Mid-Points are independent of LogRanges 
#    (Configuration C4)
ChinaldaC4 <- lda(ChinaT,ChinaTemp$GeoReg,Config=4)
cat("China linear discriminant analysis results =\n") ; print(ChinaldaC4)
cat("lda Prediction results:\n")
print(predict(ChinaldaC4,ChinaT)$class)

#Estimate error rates by three-fold cross-validation, replicated five times
CVldaC4 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,Config=4,kfold=3,CVrep=5)
summary(CVldaC4[,,"Clerr"])
glberrors <- 
	apply(CVldaC4[,,"Nk"]*CVldaC4[,,"Clerr"],1,sum)/apply(CVldaC4[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

}
