\name{Idtmclust-methods}
\docType{methods}
\alias{Idtmclust-methods}
\alias{Idtmclust,IData-method}
\alias{Idtmclust}   
\title{Methods for function Idtmclust in Package `MAINT.Data'}

\description{Performs Gaussian model based clustering for interval data}

\usage{
Idtmclust(Idt, G = 1:9, CovCase=1:4, SelCrit=c("BIC","AIC"),
  Mxt=c("Hom","Het","HomandHet"), control=EMControl()) 
}

\arguments{

  \item{Idt}{An IData object representing interval-valued entities.}
  \item{G}{An integer vector specifying the numbers of mixture components (clusters) for which the BIC is to be calculated.}
  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}
   \item{SelCrit}{The model selection criterion.}
  \item{control}{A list of control parameters for EM. The defaults are set by the call \code{EMControl()}.}
  \item{Mxt}{The type of Gaussian mixuture assumed by Idtmclust. Alternatives are \dQuote{Hom} (default) for homoceastic mixtures,
    \dQuote{Het} for heterocedastic mixtures, and \dQuote{HomandHet} for both homoceastic and heterocedastic mixtures.
  }
}

\value{
 An object of class \code{\linkS4class{IdtMclust}} providing the optimal (according to BIC) mixture model estimation.
}

\section{References}{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Brito, P., Duarte Silva, A. P. and Dias, J. G. (2015), Probabilistic Clustering of Interval Data. \emph{Intelligent Data Analysis} \bold{19}(2), 293--313.\cr

Fraley, C., Raftery, A. E., Murphy, T. B. and Scrucca, L.  (2012), mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation. 
Technical Report No. 597, Department of Statistics, University of Washington.  
}

\seealso{
  \code{\linkS4class{IdtMclust}},  \code{\link{EMControl}}, \code{\linkS4class{EMControl}} 
}

\examples{

# Create an Interval-Data object containing the intervals of loan data 
# (from the Kaggle Data Science platform) aggregated by loan purpose

LbyPIdt <- IData(LoansbyPurpose_minmaxDt,
  VarNames=c("ln-inc","ln-revolbal","open-acc","total-acc")) 

print(LbyPIdt)

\dontrun{

#Fit homocedastic Gaussian mixtures with up to four components

mclustres <- Idtmclust(LbyPIdt,G=1:4)
print(mclustres)

#Display the results of the best mixture according to the BIC

summary(mclustres,parameters=TRUE,classification=TRUE)
}

}

\keyword{interval data}
\keyword{model based clustering}

