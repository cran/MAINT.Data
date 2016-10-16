msnCP.dev <- function(param, y, grpind, Config, n=ifelse(is.matrix(y),nrow(y),1),
                     p=ifelse(is.matrix(y),ncol(y),length(y)), k=1,		     	
                     trace=FALSE, c2tol=1e-3, ldRtol=log(1e-5)+1-p,
                     PenF=1e12, PenC=0., nopenalty=FALSE)
{
  if (!all(is.finite(param))) return(.Machine$double.xmax)
  res <- .Call( "msnCP_dev", param, y, grpind, Config, n, p, k, trace,  
    c2tol, ldRtol, PenF, PenC, nopenalty, .Machine$double.eps, PACKAGE = "MAINT.Data" )
  if (!is.finite(res))  return(.Machine$double.xmax/10)
  res
}

