ConfMat <- function(origcl, predcl, otp=c("absandrel","abs","rel"), dec=3)
{
  otp <- match.arg(otp)
  
  if (!is.factor(origcl)) stop("Argument",origcl,"is not a factor\n")
  if (!is.factor(predcl)) stop("Argument",predcl,"is not a factor\n")
  if (length(origcl)!=length(predcl)) 
    stop("Lengths of arguments",origcl,"and",predcl,"do not agree with each other\n")
  clvls <- levels(origcl)
  if (any(clvls != levels(predcl))) stop("Levels of arguments",origcl,"and",predcl,"do not agree with each other\n")
  k <- length(clvls)

  Cmat <- matrix(nrow=k,ncol=k,dimnames=list(clvls,clvls))
  
  for (i in 1:k) for (j in 1:k) Cmat[i,j] <- length(which(origcl==clvls[i] & predcl==clvls[j]))
  if (otp!="abs") Cmat1 <- Cmat/rowSums(Cmat) 
  
  if (otp=="absandrel") {
    cat("Absolute frequencies (original classes in rows, predicted in columns)\n\n")
    print(Cmat)
    cat("\nRelative frequencies (original classes in rows, predicted in columns)\n\n")
    print(round(Cmat1,dec))
    invisible(list(AbsFreq=Cmat,RelFreq=Cmat1))
  }  
  else if (otp=="abs") return(Cmat)  
  else if (otp=="rel") {
    cat("Relative frequencies (original classes in rows, predicted in columns)\n")
    print(round(Cmat1,dec))
    invisible(Cmat1)
  }  
}
