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

\S4method{plot}{IData}(x, y, \dots)
\S4method{plot}{IData}(x, \dots)

}

\arguments{
  \item{x}{An object of type \linkS4class{IData} representing the values of an Interval-value variable.}

  \item{y}{An object of type \linkS4class{IData} representing the values of a second Interval-value variable, to be displayed along y (vertical) coordinates.}

  \item{\dots}{ Graphical arguments to be passed to methods. Currently the following arguments are allowed: }
  \describe{
    \item{type}{What type of plot should de drawn. Alternatives are "crosses" (default) and "rectangle".}
    \item{main}{An overall title for the plot. }
    \item{sub}{A subtitle for the plot. }
    \item{xlab}{A title for the x axis.}
    \item{ylab}{A title for the y axis.}
  }
}

\seealso{
  \linkS4class{IData} 
}

