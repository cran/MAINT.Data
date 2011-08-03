\name{MAINT.Data-package}
\alias{MAINT.Data-package}
\alias{MAINT.Data}
\docType{package}
\title{Modelling and Analizing Interval Data}
\description{MAINT.DATA implements methodologies for modelling Interval Data, considering five different possible configurations
structures for the variance-covariance matrix.
It introduces a data class for representing interval data and
includes functions and methods for parametric modelling and analysing of interval data.
It performs maximum likelihood estimation and statistical tests as well as (M)ANOVA and Linear and Quadratic Discriminant Analysis for all considered configurations.
}
\details{In the classical model of multivariate data analysis, data is represented in a data-array where n ``individuals" (usually in rows) take exactly
one value for each variable (usually in columns).
Symbolic Data Analysis (see, e.g., Noirhomme-Fraiture and Brito (2011)) provides a framework where new variable types allow to take directly into account variability and/or uncertainty associated to each single ``individual",
by allowing multiple, possibly weighted, values for each variable. 
New variable types - interval, categorical multi-valued and modal variables - have been introduced.\cr
We focus on the analysis of interval data, i.e., where
elements are described by variables whose values are intervals.
Parametric inference methodologies based on probabilistic models for interval variables are developed in Brito and Duarte Silva (2011) where each interval is represented by its midpoint and log-range,for which Normal and Skew-Normal distributions are assumed.
The intrinsic nature of the interval variables leads to special structures of the variance-covariance matrix, which are represented by five different possible configurations.\cr
MAINT.DATA implements the proposed methodologies in R,
introducing a data class for representing interval data; it
includes functions for modelling and analysing interval data, 
in particular maximum likelihood estimation and statistical tests for the different considered configurations. 
Methods for (M)ANOVA and Linear and Quadratic Discriminant Analysis of this data class are also provided.

\tabular{ll}{
Package: \tab MAINT.Data\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2011-08-03\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
}
\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>

Maintainer: Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
}

\references{
Brito, P., Duarte Silva, A. P.  (2011): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics (forthcoming).\cr

Noirhomme-Fraiture, M., Brito, P. (2011): "Far Beyond the Classical Data Models: Symbolic Data Analysis". Statistical Analysis and Data Mining, Volume 4, Issue 2, 157-170.
}

\keyword{MAINT.Data}
\keyword{Interval Data}
\keyword{Symbolic Data Analysis}
\keyword{MANOVA for Interval Data}
\keyword{Parametric modelling of Interval Data}
\keyword{Statistical tests for Interval Data}

\examples{
# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#Display the first and last observations
head(ChinaT)
tail(ChinaT)

#Print summary statistics
summary(ChinaT)

# Estimate parameters by maximum likelihood, assuming 
# the classical (unrestricted) covariance configuration C1
ChinaEC1 <- mle(ChinaT,Config=1)
cat("China maximum likelhiood estimation results:\n") ; print(ChinaEC1)
cat("Standard Errors of Estimators:\n") ; print(stdEr(ChinaEC1))

# Estimate parameters by maximum likelihood, 
# assuming that one of the restricted covariance configurations holds
ChinaEC2345 <- mle(ChinaT,Config=2:5)
cat("China maximum likelhiood estimation results:\n") ; print(ChinaEC2345)
cat("Standard Errors of Estimators:\n") ; print(stdEr(ChinaEC2345))

#MANOVA tests  
ManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg)
cat("China, MANOVA by geografical regions results:\n") ; print(ManvChina)

#Linear Discriminant Analysis
Chinalda <- lda(ManvChina)
cat("China linear discriminant analysis results =\n") ; print(Chinalda)
cat("lda Prediction results:\n")
print(predict(Chinalda,ChinaT)$class)

#Estimate error rates by ten-fold cross-validation
CVlda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,
Config=BestModel(ManvChina@H1res),CVrep=1)
summary(CVlda[,,"Clerr"])
glberrors <- 
  apply(CVlda[,,"Nk"]*CVlda[,,"Clerr"],1,sum)/apply(CVlda[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

}
