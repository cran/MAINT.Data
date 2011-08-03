\name{BestModel-methods}
\docType{methods}
\alias{BestModel-methods}
\alias{BestModel}   
\title{ Methods for Function BestModel in Package `MAINT.Data'}

\description{Selects the best model according to the chosen selection criterion (currently, AIC or BIC) }

\usage{BestModel(Idt,SelCrit=c("IdtCrt","AIC","BIC"))}

\arguments{
  \item{Idt}{An IData object representing NObs interval-valued entities.}

  \item{SelCrit}{The model selection criterion}
}

\value{An integer with the index of the model chosen by the selection criterion}

\keyword{methods}
\keyword{AIC}
\keyword{BIC}
\keyword{interval data}


