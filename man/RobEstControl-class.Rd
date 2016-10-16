\name{RobEstControl-class}
\docType{class}
\alias{RobEstControl-class}

\title{Class 'RobEstControl' - contains control parameters for the robust estimation of parametric interval data models.}
\description{
    This class extends the \code{CovControlMcd} class
    and contains control parameters for the robust estimation of parametric interval data models.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RobEstControl", ...)}
or by calling the constructor-function \code{RobEstControl}.
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Inherited from class \code{"CovControlMcd"}. Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*Idt@NIVar observations are used for computing the trimmed likelihood. Allowed values are between 0.5 and 1.}
    \item{\code{nsamp}:}{Inherited from class \code{"CovControlMcd"}. Number of subsets used for initial estimates. Note that when argument \sQuote{getalpha} is set to \dQuote{TwoStep} the final value of \sQuote{alpha} is estimated by a two-step procedure and the value of argument \sQuote{alpha} is only used to specify the size of the samples used in the first step.}
    \item{\code{scalefn}:}{Inherited from class \code{"CovControlMcd"} and not used in the package \\sQuote{Maint.Data.}}
    \item{\code{maxcsteps}:}{Inherited from class \code{"CovControlMcd"} and not used in the package \\sQuote{Maint.Data.}}
    \item{\code{seed}:}{Inherited from class \code{"CovControlMcd"}. Starting value for random generator. Default is \code{seed = NULL.}}
    \item{\code{use.correction}:}{Inherited from class \code{"CovControlMcd"}. Whether to use finite sample correction factors.
Default is \code{use.correction=TRUE}.}
    \item{\code{trace}, \code{tolSolve}:}{Inherited from class \code{"CovControl"}.} 
    \item{\code{ncsteps}:}{The maximum number of concentration steps used each iteration of the fasttle algorithm.}
    \item{\code{getalpha}:}{Argument specifying if the \sQuote{alpha} parameter (roughly the percentage of the sample used for computing the trimmed likelihood) should be estimadted from the data, or if the value of the argument \sQuote{alpha} should be used instead. When set to \dQuote{TwoStep}, \sQuote{alpha} is estimated by a two-step procedure with the value of argument \sQuote{alpha} specifying the size of the samples used in the first step. Otherwise, with the value of argument \sQuote{alpha} is used directly.}   
    \item{\code{getkdblstar}:}{Argument specifying the size of the initial small (in order to minimize the probability of outliers) subsets. If set to the string \dQuote{Twopplusone} (default) the initial sets have twice the number of interval-value variables plus one (i.e., they are the smaller samples that lead to a non-singular covaraince estimate). Otherwise, an integer with the size of the initial sets.}
    \item{\code{outlin}:}{The type of outliers to be consideres. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}
    \item{\code{trialmethod}:}{The method to find a trial subset used to initialize each replication of the fasttle algorithm. The current options are \dQuote{simple} (default) that simply selects \sQuote{kdblstar} observations at random, and \dQuote{Poolm} that divides the original sample into \sQuote{m} non-overlaping subsets, applies the \sQuote{simple trial} and the refinement methods to each one of them, and merges the results into a trial subset.}
    \item{\code{m}:}{Number of non-overlaping subsets used by the trial method when the argument of \sQuote{trialmethod} is set to 'Poolm'.}
    \item{\code{reweighted}:}{Should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}
    \item{\code{otpType}:}{The amount of output returned by fasttle. Current options are \dQuote{OnlyEst} (default) where only an \sQuote{IdtE} object with the fasttle estimates is returned, \dQuote{SetMD2andEst} which returns a list with an \sQuote{IdtE} object of fasttle estimates, a vector with the final trimmed subset elements used to compute these estimates and the corresponding robust squared Mahalanobis distances, and \dQuote{SetMD2EstandPrfSt} wich returns a list with the previous three components plust a list of some performance statistics concerning the algorithm execution.}
  }
}

\section{Extends}{
    Class \code{"\linkS4class{CovControlMcd}"}, directly.
    Class \code{"\linkS4class{CovControl}"} by CovControlMcd, distance 2.
}

\section{Methods}{
No methods defined with class "RobEstControl" in the signature.
}

\seealso{
  \code{\link{RobEstControl}}, \code{\link{fasttle}}, \code{\link{RobMxtDEst}}, \code{\link{Roblda}}, \code{\link{Robqda}}
}

\references{
Hadi, A. S. and Luceno, A. (1997), Maximum trimmed likelihood estimators: a unified approach, examples, and algorithms.
\emph{Computational Statistics and Data Analysis} \bold{25}(3), 251--272.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}


\keyword{classes}
\keyword{robust}
\keyword{interval data}
