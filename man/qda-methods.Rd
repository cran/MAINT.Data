\name{qda-methods}
\docType{methods}
\alias{qda}
\alias{qda-methods}
\alias{qda,IData-method}
\alias{qda,IdtHetNMANOVA-method}
\title{Methods for Function qda in Package `MAINT.Data'}
\description{Quadratic Discriminant Analysis}

\section{Methods}{
\describe{
\item{\code{signature(Idt = "IData")}}{Performs Quadratic Discriminant Analysis on interval-valued data}
}
}

\keyword{methods}
\keyword{Quadratic Discrimiant Analysis}
\keyword{interval data}

\usage{

\S4method{qda}{IdtHetNMANOVA}(x,prior="proportions",selmodel=BestModel(H1res(x)),\dots)

\S4method{qda}{IData}(x, grouping, prior="proportions", tol=1.0e-4, 
subset=1:nrow(x) ,Config=1:5, SelCrit=c("AIC","BIC"))
}

\arguments{

  \item{x}{An IdtHetNMANOVA object representing the results of a (heterocedastic normal-based) MANOVA of Interval Data, or an IData object representing NObs interval-valued entities}

  \item{grouping}{Factor indicating the group to which each observation belongs to}

  \item{prior}{The prior probabilities of class membership.  If unspecified,
          the class proportions for the training set are used.  If
          present, the probabilities should be specified in the order
          of the factor levels}

  \item{selmodel}{The index of model of IdtClMANOVA object to be used in the construction of the discriminant functions}

  \item{tol}{A tolerance to decide if a matrix is singular; it will reject
          variables and linear combinations of unit-variance variables
          whose standard deviation is less than tol}

  \item{subset}{An index vector specifying the cases to be used in the
          training sample.  (NOTE: If given, this argument must be
          named}

   \item{Config}{Configuration of the variance-covariance matrix: a vector of integers between 1 and 5}

   \item{SelCrit}{The model selection criterion}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{An object of type Idtqda, containing the quadratic discriminant functions}

\seealso{ \code{\linkS4class{Idtlda}}, \code{\link{DACrossVal}}, \code{\link{lda}} }

\examples{

#Create an Interval-Data object containing the intervals of temperatures 
# by quarter for 60 Chinese meteorological stations.
ChinaT <- IData(ChinaTemp[1:8])

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.
ChinaT <- IData(ChinaTemp[1:8])

#Quadratic Discriminant Analysis, assuming independent Interval Variables 
#   (Configuration C3)
ChinaqdaC3 <- qda(ChinaT,ChinaTemp$GeoReg,Config=3)
cat("China quadratic discriminant analysis results =\n") ; print(ChinaqdaC3)
cat("qda Prediction results:\n")
print(predict(ChinaqdaC3,ChinaT)$class)

#Estimate error rates by three-fold cross-validation, replicated five times
CVqdaC3 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=qda,Config=3,kfold=3,CVrep=5)
summary(CVqdaC3[,,"Clerr"])
glberrors <- 
	apply(CVqdaC3[,,"Nk"]*CVqdaC3[,,"Clerr"],1,sum)/apply(CVqdaC3[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

}
