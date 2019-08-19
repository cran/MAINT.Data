\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,IData-method}

\title{Methods for function plot in Package `MAINT.Data'}
\description{S4 methods for function plot. As in the generic plot S3 `graphics' method, these methods plot Interval-valued data contained in \linkS4class{IData} objects.}

\keyword{methods}
\keyword{plot}
\keyword{Interval Data}

\usage{

\S4method{plot}{IData}(x, y, type=c("crosses","rectangles"), append=FALSE, \dots)
\S4method{plot}{IData}(x, casen=NULL, layout=c("vertical","horizontal"), append=FALSE, \dots)

}

\arguments{
  \item{x}{An object of type \linkS4class{IData} representing the values of an Interval-value variable.}

  \item{y}{An object of type \linkS4class{IData} representing the values of a second Interval-value variable, to be displayed along y (vertical) coordinates.}

    \item{type}{What type of plot should de drawn. Alternatives are "crosses" (default) and "rectangles".}

    \item{append}{A boolean flag indicating if the interval-valued variables should be displayed in a new plot, or added to an existing plot.}

    \item{casen}{An optional character string with the case names.}

    \item{layout}{The axes along which the interval-valued variables be displayed. Alternatives are "vertical" (default) and "horizontal".}
    
  \item{\dots}{ Graphical arguments to be passed to methods.}
}

\seealso{
  \linkS4class{IData} 
}


