\name{DACrossVal}

\title{Cross Validation for Discriminant Analysis Classification Rules}

\alias{DACrossVal}

\description{
\sQuote{DACrossVal} evaluates the performance of a Discriminant Analysis training sample algorithm by k-fold 
Cross-Validation.  
}

\usage{DACrossVal(data, grouping, TrainAlg, EvalAlg=EvalClrule, 
Strfolds=TRUE, kfold=10, CVrep=20, prior="proportions", \dots)}

\arguments{

  \item{data}{Matrix, data frame or Interval Data object of observations.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{TrainAlg}{A function with the training algorithm. It should return an object that can be used as input to the argument of \sQuote{EValAlg}.}

  \item{EvalAlg}{A function with the evaluation algorithm. By default set to \sQuote{EvalClrule} which returns a list with components \dQuote{err} (estimates of error rates by class) and \dQuote{Nk} (number of out-sample observations by class). This default can be used for all \sQuote{TrainAlg} 
arguments that return an object with a predict method returning a list with a \sQuote{class} component (a factor) containing the classification results.}

  \item{Strfolds}{Boolean flag indicating if the folds should be stratified according to the original class
proportions (default), or randomly generated from the whole training sample, ignoring class membership.}

  \item{kfold}{Number of training sample folds to be created in each replication.}

  \item{CVrep}{Number of replications to be performed.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If present, the probabilities should be specified in the order of the factor levels.}

  \item{\dots}{Further arguments to be passed to \sQuote{TrainAlg} and \sQuote{EvalAlg}.} 
}

\value{A three dimensional array with the number of tested observations, and estimated classification errors for each combination of fold and replication tried. The array dimensions are defined as follows:\cr
The first dimension runs through the different fold-replication combinations.\cr
The second dimension represents the classes. \cr
The third dimension has two named levels representing respectively the number of observations tested (\dQuote{Nk}), and the estimated classification errors (\dQuote{Clerr}).}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{lda}}, \code{\link{qda}}, \code{\link{IData}}}

\examples{

# Compare performance of linear and quadratic discriminant analysis with 
#  Configurations C1 and c4 on the ChinaT data set by 5-fold cross-validation 
#  replicated twice

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Classical (configuration 1) Linear Discriminant Analysis 

CVldaC1 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,Config=1,kfold=5,CVrep=2)
summary(CVldaC1[,,"Clerr"])
glberrors <- 
	apply(CVldaC1[,,"Nk"]*CVldaC1[,,"Clerr"],1,sum)/apply(CVldaC1[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

# Linear Discriminant Analysis with configuration 4

CVldaC4 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,Config=4,kfold=5,CVrep=2)
summary(CVldaC4[,,"Clerr"])
glberrors <- 
	apply(CVldaC4[,,"Nk"]*CVldaC4[,,"Clerr"],1,sum)/apply(CVldaC4[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

# Classical (configuration 1) Quadratic Discriminant Analysis 

CVqdaC1 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=qda,Config=1,kfold=5,CVrep=2)
summary(CVqdaC1[,,"Clerr"])
glberrors <- 
	apply(CVqdaC1[,,"Nk"]*CVqdaC1[,,"Clerr"],1,sum)/apply(CVqdaC1[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

# Quadratic Discriminant Analysis with configuration 4

CVqdaC4 <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=qda,Config=4,kfold=5,CVrep=2)
summary(CVqdaC4[,,"Clerr"])
glberrors <- 
	apply(CVqdaC4[,,"Nk"]*CVqdaC4[,,"Clerr"],1,sum)/apply(CVqdaC4[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

}

\keyword{discriminant analysis}
\keyword{cross-validation}

