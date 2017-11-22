\name{EMControl}
\alias{EMControl}

\title{ Constructor function for objects of class EMControl}
\description{
  This function will create a control object of class \code{\linkS4class{EMControl}} containing the control parameters 
  for the EM algorithm  used in estimation of Gaussian mixtures by function \code{\link{Idtmclust}}.  
}

\usage{
  EMControl(nrep=100, maxiter=1000, convtol=0.01, protol=1e-6, seed=NULL)
}

\arguments{
  \item{nrep}{Number of replications (different randomly generated starting points) of the EM algorithm.}
  \item{maxiter}{Maximum number of iterations in each replication of the EM algorithm.}
  \item{convtol}{Numeric tolerance for testing the convergence of the EM algorithm. Convergence is assumed when the log-likelihood changes less than convtol.}
  \item{protol}{Numeric tolerance for the mixture proportions. Proportions below protol, considered to be zero, are not allowed.}
  \item{seed}{Starting value for random generator.}
}

\value{An \code{EMControl} object}

\seealso{
  \code{\link{Idtmclust}}
}

