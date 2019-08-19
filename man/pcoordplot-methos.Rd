\name{pcoordplot-methods}
\docType{methods}
\alias{pcoordplot}
\alias{pcoordplot-methods}
\alias{pcoordplot,IdtMclust-method}

\title{Parallel coordinates plot.}
\description{Method pcoordplot displays a parallel coordinates plot, representing the results stored in an IdtMclust-method object.}

\keyword{methods}
\keyword{Parallel coordinates plot}
\keyword{Interval Data}
\keyword{model based clustering}

\usage{

\S4method{pcoordplot}{IdtMclust}(x,title="Parallel Coordinate Plot",
Seq=c("AllMidP_AllLogR","MidPLogR_VarbyVar"), G=BestG(x), \dots)

}

\arguments{
  \item{x}{An object of type \dQuote{IdtMclust} representing the the clusterig results of an Interval-valued data set obtainde by the function \dQuote{IdtMclust}.}

  \item{title}{The title of the plot.}

  \item{\dots}{ Graphical arguments to be passed to methods}

  \item{Seq}{The ordering of the coordinates in the plot. Available options are:\cr
\dQuote{AllMidP_AllLogR}: all MidPoints followed all LogRanges, in the same variable order.\cr
\dQuote{MidPLogR_VarbyVar}: MidPoints followed by LogRanges, variable by variable.\cr
}

  \item{G}{An integer vector specifying the numbers of mixture components (clusters) of the solution to be displayed.}

}

\seealso{
  \linkS4class{IdtMclust}, \code{\link{Idtmclust}} 
}


