\name{da-methods}
\docType{methods}
\alias{lda-methods}
\alias{lda,IData-method}
\alias{lda,IdtMxtNDE-method}
\alias{lda,IdtClMANOVA-method}
\alias{lda,IdtLocNSNMANOVA-method}
\alias{qda-methods}
\alias{qda,IData-method}
\alias{qda,IdtMxtNDE-method}
\alias{qda,IdtHetNMANOVA-method}
\alias{qda,IdtGenNSNMANOVA-method}
\alias{Roblda}
\alias{Roblda-methods}
\alias{Roblda,IData-method}
\alias{Robqda}
\alias{Robqda-methods}
\alias{Robqda,IData-method}
\alias{snda-methods}
\alias{snda,IData-method}
\alias{snda,IdtLocSNMANOVA-method}
\alias{snda,IdtLocNSNMANOVA-method}
\alias{snda,IdtGenSNMANOVA-method}
\alias{snda,IdtGenNSNMANOVA-method}

\title{Discriminant Analysis of Interval Data}

\description{lda and qda perform linear and quadratic discriminant analysis of Interval Data based on classic estimates of a mixture of  Gaussian models. Roblda and Robqda do the same using robust estimates of location and scatter. snda performs discriminant analysis of Interval Data based on estimates of mixtures of Skew-Normal models}

\usage{
\S4method{lda}{IData}( x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10, 
  subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, \dots )

\S4method{lda}{IdtMxtNDE}(x, prior="proportions", selmodel=BestModel(x), egvtol=1.0e-10,
  silent=FALSE, \dots )

\S4method{lda}{IdtClMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  egvtol=1.0e-10, silent=FALSE, \dots )

\S4method{lda}{IdtClMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  egvtol=1.0e-10, silent=FALSE, \dots )

\S4method{lda}{IdtLocNSNMANOVA}( x, prior="proportions", 
  selmodel=BestModel(H1res(x)@NMod), egvtol=1.0e-10, silent=FALSE, \dots )

\S4method{qda}{IData}( x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
  CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, \dots )

\S4method{qda}{IdtMxtNDE}(x, prior="proportions", selmodel=BestModel(x), silent=FALSE, 
  \dots )

\S4method{qda}{IdtHetNMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)), 
  silent=FALSE, \dots )

\S4method{qda}{IdtGenNSNMANOVA}( x, prior="proportions", 
  selmodel=BestModel(H1res(x)@NMod), silent=FALSE, \dots )

\S4method{Roblda}{IData}( x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10,
  subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, 
  CovEstMet=c("Pooled","Globdev"), SngDMet=c("fasttle","fulltle"),
  Robcontrol=RobEstControl(), \dots )

\S4method{Robqda}{IData}( x, grouping, prior="proportions", CVtol=1.0e-5, 
  subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE,
  SngDMet=c("fasttle","fulltle"), Robcontrol=RobEstControl(), \dots )

\S4method{snda}{IData}(x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
  CovCase=1:4, SelCrit=c("BIC","AIC"), Mxt=c("Loc","Gen"), \dots )

\S4method{snda}{IdtLocSNMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  egvtol=1.0e-10, silent=FALSE, \dots )

\S4method{snda}{IdtLocNSNMANOVA}( x, prior="proportions",
  selmodel=BestModel(H1res(x)@SNMod), egvtol=1.0e-10, silent=FALSE, \dots )

\S4method{snda}{IdtGenSNMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  silent=FALSE, \dots )

\S4method{snda}{IdtGenNSNMANOVA}( x, prior="proportions",
  selmodel=BestModel(H1res(x)@SNMod), silent=FALSE, \dots )

}

\arguments{

  \item{x}{An object of class \code{\linkS4class{IData}}, \code{\linkS4class{IdtLocSNMANOVA}}, \code{\linkS4class{IdtLocNSNMANOVA}},\code{\linkS4class{IdtGenSNMANOVA}} or \code{\linkS4class{IdtGenNSNMANOVA}} with either the original Interval Data, or the results of a Interval Data Skew-Normal MANOVA, from which the discriminant analysis will be based.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

  \item{subset}{An index vector specifying the cases to be used in the analysis.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

  \item{SelCrit}{The model selection criterion.}

  \item{silent}{A boolean flag indicating wether a warning message should be printed if the method fails.}
 
  \item{CovEstMet}{Method used to estimate the common covariance matrix in \code{Roblda} (Robust linear discriminant analysis). Alternatives are \dQuote{Pooled} (default) for a pooled average of the the robust within-groups covariance estimates, and \dQuote{Globdev} for a global estimate based on all deviations from the groups multivariate l1 medians. See Todorov and Filzmoser (2009) and \code{pcaPP.l1median} for details.}
 
  \item{SngDMet}{Algorithm used to find the robust estimates of location and scatter. Alternatives are \dQuote{fasttle} (default) and \dQuote{fulltle}.}
 
  \item{Robcontrol}{A control object (S4) of class \code{\link{RobEstControl-class}}
    containing estimation options - same as these provided in the function
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
 
  \item{Mxt}{Indicates the type of mixing distributions to be considered. Current alternatives are \dQuote{Hom} (homocedastic) and \dQuote{Het} (hetereocedasic) for Gaussian models, \dQuote{Loc} (location model -- groups differ only on their location parameters) and \dQuote{Gen} \dQuote{Loc} (general model -- groups differ on all parameters) for Skew-Normal models.}

  \item{selmodel}{Selected model from a list of candidate models saved in object x.}

  \item{egvtol}{Tolerance level for the eigenvalues of the product of the inverse within by the between covariance matrices. When a  eigenvalue has an absolute value below egvtol, it is considered to be zero.}

  \item{\dots}{Other named arguments.}
}

\references{
Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based   approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}

\seealso{\code{\linkS4class{IData}}, \code{\linkS4class{IdtLocSNMANOVA}},\code{\linkS4class{IdtLocSNMANOVA}},\code{\linkS4class{IdtLocSNMANOVA}},\code{\linkS4class{IdtLocSNMANOVA}}, \code{pcaPP.l1median}.}

\keyword{methods}
\keyword{Interval Data}
\keyword{Discriminant Analysis Interval Data}

\examples{
# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#Linear Discriminant Analysis

ChinaT.lda <- lda(ChinaT,ChinaTemp$GeoReg)
cat("Temperatures of China -- linear discriminant analysis results:\n")
print(ChinaT.lda)
cat("lda Prediction results:\n")
print(predict(ChinaT.lda,ChinaT)$class)

\dontrun{
#Estimate error rates by ten-fold cross-validation replicated 20 times  

CVlda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,CovCase=CovCase(ChinaT.lda))
summary(CVlda[,,"Clerr"])
glberrors <- 
  apply(CVlda[,,"Nk"]*CVlda[,,"Clerr"],1,sum)/apply(CVlda[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")
}

#Quadratic Discriminant Analysis

ChinaT.qda <- qda(ChinaT,ChinaTemp$GeoReg)
cat("Temperatures of China -- qda discriminant analysis results:\n")
print(ChinaT.qda)

\dontrun{
#Estimate error rates by ten-fold cross-validation replicated 20 times  

CVqda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=qda,CovCase=CovCase(ChinaT.qda))
summary(CVqda[,,"Clerr"])
glberrors <- 
  apply(CVqda[,,"Nk"]*CVqda[,,"Clerr"],1,sum)/apply(CVqda[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

# Skew-Normal based discriminant analysis, asssuming that the different regions may differ
# in all SkewNormal parameters

cat("Temperatures of China -- SkewNormal general model discriminant analysis results:\n")
ChinaT.gensnda <- snda(ChinaT,ChinaTemp$GeoReg,Mxt="Gen")
print(ChinaT.gensnda)

#Estimate error rates by three-fold cross-validation without replication  

CVgensnda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=snda,Mxt="Gen",
  CovCase=CovCase(ChinaT.gensnda),kfold=3,CVrep=1)
summary(CVgensnda[,,"Clerr"])
glberrors <- 
  apply(CVgensnda[,,"Nk"]*CVgensnda[,,"Clerr"],1,sum)/apply(CVgensnda[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

#Robust Quadratic Discriminant Analysis

ChinaT.rqda <- Robqda(ChinaT,ChinaTemp$GeoReg)
cat("Temperatures of China -- robust qda discriminant analysis results:\n")
print(ChinaT.rqda)

#Estimate error rates by ten-fold cross-validation with 5 replications 

CVrqda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=Robqda,CovCase=CovCase(ChinaT.rqda),
   CVrep=5)
summary(CVrqda[,,"Clerr"])
glberrors <- 
  apply(CVrqda[,,"Nk"]*CVrqda[,,"Clerr"],1,sum)/apply(CVrqda[,,"Nk"],1,sum)
cat("Average global classification error =",mean(glberrors),"\n")

}

}
