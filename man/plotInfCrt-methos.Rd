\name{plotInfCrt-methods}
\docType{methods}
\alias{plotInfCrt}
\alias{plotInfCrt-methods}
\alias{plotInfCrt,IdtMclust-method}

\title{Information criteria plot.}
\description{Method plotInfCrt displays a plot representing the values of an approprtiate information criterion (currently either BIC or AIC) for the models whose results are stored in an IdtMclust-method object. A supplementary short output message prints the values of the chosen criterion for the 'nprin' best models.}

\keyword{methods}
\keyword{BIC}
\keyword{AIC}
\keyword{Interval Data}
\keyword{model based clustering}

\usage{

\S4method{plotInfCrt}{IdtMclust}(object,crt=object@SelCrit,legpos="bottomleft",nprnt=5, \dots)

}

\arguments{
  \item{object}{An object of type \dQuote{IdtMclust} representing the the clusterig results of an Interval-valued data set obtainde by the function \dQuote{IdtMclust}.}

  \item{crt}{The information criteria whose values are to be displayed.}

  \item{legpos}{Legend position. Alternatives are \dQuote{bottomleft} (default), \dQuote{bottomright}, \dQuote{topleft} and \dQuote{bottomright}.} 

  \item{nprnt}{Number of solutions for which the value of the information criterio should be printed in an suplmentary short output message.}

  \item{\dots}{ Graphical arguments to be passed to methods.}

}

\seealso{
  \linkS4class{IdtMclust}, \code{\link{Idtmclust}} 
}


