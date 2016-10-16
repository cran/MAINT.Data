msnCP.dev.grad <- function(param, y, grpind, Config, n=ifelse(is.matrix(y),nrow(y),1),
  p=ifelse(is.matrix(y),ncol(y),length(y)), k=1, inoptm=TRUE, 		     	
  trace=FALSE, c2tol=1e-3, ldRtol=log(1e-5)+1-p, beta0tol=1e-6, PenF=1e12, ...)
{
  if (!all(is.finite(param))) {
    if (inoptm)  {
      return(rep(0.,length(param)))
    } else {
      return(NULL)
    }
  }
  if (!is.matrix(y)) { dim(y) <- c(1,p) }

  grad <- .Call( "msnCP_dev_grad", param, y, grpind, Config, n, p, k, 
    trace, c2tol, ldRtol, beta0tol, PenF, .Machine$double.eps, PACKAGE = "MAINT.Data" )

  if (!all(is.finite(grad))) { 
    if (inoptm)  {
      return(rep(0.,length(param)))
    } else {
      return(NULL)
    }
  }

  grad
}

