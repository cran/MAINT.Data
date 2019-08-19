CheckMSing <- function(M,limlnk2,limdiag=NULL)
{
   if (!is.null(limdiag) && any(diag(M)<limdiag)) return(TRUE)
   trM <- sum(diag(M)) 
   lnMdet <- determinant(M,logarithm=TRUE)
   if (lnMdet$sign==-1) return(TRUE) 
   p <- ncol(M)
   if ( p*log(trM) - lnMdet$modulus < limlnk2) return(FALSE)   
   Megval <- eigen(M,only.values=TRUE)$values
   if (Megval[p] < 0.) return(TRUE)
   if (log(Megval[1]/Megval[p]) < limlnk2) return(FALSE)
   else return(TRUE)
}

CheckSigmaSing <- function(Config,Sigma,limlnk2=log(1e15),limdiag=1e-300)
{
    if (any(diag(Sigma)<limdiag)) return(TRUE)
    if (Config==1) return(CheckMSing(Sigma,limlnk2))
    p <- ncol(Sigma)
    q <- p/2

    if (Config==3) {
      for (j in 1:q) {
        a <- Sigma[j,j]
        b <- Sigma[q+j,q+j]
        c <- Sigma[j,q+j]
        trSigj <- a + b
        detSigj <- a*b - c^2
        d <- sqrt(trSigj^2-4*detSigj)
        k2j <- (trSigj+d)/(trSigj-d)
        if (!is.finite(k2j) || k2j<limdiag) return(TRUE) 
        if (log(k2j) >  limlnk2) return(TRUE)
      }
      return(FALSE)
    }

    if (Config==4) {
      if (CheckMSing(Sigma[1:q,1:q],limlnk2)==TRUE) return(TRUE)
      return(CheckMSing(Sigma[(q+1):p,(q+1):p],limlnk2))
    }

    if (Config==5) {
      vars <- diag(Sigma) 
      egval1 <- max(vars)
      egvalp <- min(vars)
      if (egvalp < 0 || log(egval1/egvalp) >=  limlnk2) return(TRUE)
      else return(FALSE)
    }

}

CheckSigmakSing <- function(Config,Sigmak,limlnk2=log(1e15))
{
   k <- dim(Sigmak)[3]
   for (g in 1:k) 
     if (CheckSigmaSing(Config,Sigmak[,,g],limlnk2)==TRUE) return(TRUE)
   return(FALSE)
}

RepEMGauss <- function(X,n,p,k,ISol,pertub,nrep,Config,Homoc,maxiter,convtol,protol,Ifact=1e-3)
{
  genunif <-function(u0,amp,ubnd,EPS=1E-6) return(max(EPS,min(ubnd,u0+runif(1,0.,amp)-amp/2)))
  mychol <-function(Sigma,p,pert=1E-6) return(chol(Sigma+pert*diag(1.,p)))
  
  Vnames <- names(X)
  Onames <- rownames(X) 
  Lik <- array(nrep)
  if (Homoc==FALSE) Sigmak <- array(dim=c(p,p,k),dimnames=list(Vnames,Vnames,NULL))
  if (!is.null(pertub$tau)) tau <- array(dim=k)
  if (!is.null(pertub$z)) z <- array(dim=c(n,k))

  BestLnLik <- ISol$LnLik
  q <- p/2
  for (rep in 1:nrep)  {
    if (!is.null(pertub$z)) {
      tau <- muk <- Sigma <- Sigmak <- NULL
      ampl <- rep(1.,n)
      for (i in 1:(k-1)) {
        z[,i] <- sapply(z[,i],genunif,amp=pertub$z*ampl,ubnd=ampl)
        ampl <- ampl - z[,i]
      }
      z[,k] <- ampl
      LnLik <- MClusLikz(z,X,n,p,k,Config,Homoc)
    } else {
      z <- NULL
      if (!is.null(pertub$tau)) {
        ampl <- 1.
        for (i in 1:(k-1)) {
          tau[i] <- max(0.,genunif(ISol$tau[i],amp=pertub$tau[i]*ampl,ubnd=ampl))
          ampl <- max(0.,ampl-tau[i])
        }
        tau[k] <- ampl
      } else {
        tau <- ISol$tau
      }  
      if (!is.null(pertub$muk)) {
        muk <- matrix(rnorm(k*p,ISol$muk,pertub$muk),nrow=p,ncol=k,dimnames=list(Vnames,NULL))
      } else {
        muk <- ISol$muk
      } 
      
      if (!is.null(pertub$Stdev)) {
        q0 <- p*(p-1)/2
       	if (Homoc==TRUE) {
          if (CheckSigmaSing(Config,ISol$Sigma)==TRUE) {
            diag(ISol$Sigma) <- (1.+Ifact)*diag(ISol$Sigma)
          }
       	  Sk <- 1
       	} else {
          if (CheckSigmakSing(Config,ISol$Sigmak)==TRUE) {
            for (g in 1:k) diag(ISol$Sigmak[,,k]) <- (1.+Ifact)*diag(ISol$Sigmak[,,k])
          }
       	  Sk <- k
       	}  
        if (Config==3)  RSr <- vector("list",q)
        else if (Config==4) RSr <- vector("list",2) 
        for (g in 1:Sk)  {
          if (Homoc==TRUE) {
            Stdev <- sqrt(diag(ISol$Sigma))
            if (Config==1 || Config==2) RSr <- mychol(cov2cor(ISol$Sigma),p)
            else if (Config==3) for (v in 1:q) RSr[[v]] <- mychol(cov2cor(ISol$Sigma[c(v,q+v),c(v,q+v)]),2)
            else if (Config==4) { 
              RSr[[1]] <- mychol(cov2cor(ISol$Sigma[1:q,1:q]),q)
              RSr[[2]] <- mychol(cov2cor(ISol$Sigma[q+1:q,q+1:q]),q)
            }
          } else {
            Stdev <- sqrt(diag(ISol$Sigmak[,,g]))
            if (Config==1 || Config==2) RSr <- mychol(cov2cor(ISol$Sigmak[,,g]),p)
            else if (Config==3) for (v in 1:q) RSr[[v]] <- mychol(cov2cor(ISol$Sigmak[,,g][c(v,q+v),c(v,q+v)]),2)
            else if (Config==4) { 
              RSr[[1]] <- mychol(cov2cor(ISol$Sigmak[,,g][1:q,1:q]),q)
              RSr[[2]] <- mychol(cov2cor(ISol$Sigmak[,,g][q+1:q,q+1:q]),q)
            }
          } 
          Stdev <- abs(rnorm(p,Stdev,pertub$Stdev[(g-1)*p+1:p]))
          if (Config==5) {
            if (Homoc==TRUE) Sigma <- diag(Stdev^2)
            else Sigmak[,,g] <- diag(Stdev^2)
          } else {
            if (Config==1 || Config==2) {
              uptrigel <- upper.tri(RSr) 
              RSr[uptrigel] <- rnorm(q0,RSr[uptrigel],pertub$cor[(g-1)*q0+1:q0])
              R <- cov2cor(t(RSr) %*% RSr)        
              if (Homoc==TRUE) {
                Sigma <- R * outer(Stdev,Stdev)
              } else {
                Sigmak[,,g] <- R * outer(Stdev,Stdev)
              }
            }  else if (Config==3 || Config==4) {            
              if (Homoc==TRUE) Sigma <- matrix(0.,nrow=p,ncol=p)
              else Sigmak[,,g] <- 0.
              if (Config==3) {            
                for (v in 1:q) {            
                  RSr[[v]][1,2] <- rnorm(1,RSr[[v]][1,2],pertub$cor[(g-1)*q+v])
                  R <- cov2cor(t(RSr[[v]]) %*% RSr[[v]])
                  indices <- c(v,q+v)         
                  if (Homoc==TRUE) {
                    Sigma[indices,indices] <- R * outer(Stdev[indices],Stdev[indices])
                  } else {
                    Sigmak[indices,indices,g] <- R * outer(Stdev[indices],Stdev[indices])
                  }
                }
              } else if (Config==4) {     
                q4 <- q*(q-1)/2       
                uptrigel <- upper.tri(RSr[[1]])  # it could just as well be upper.tri(RSr[[2]])... 
                for (b in 1:2) {            
                  RSr[[b]][uptrigel] <- rnorm(q4,RSr[[b]][uptrigel],pertub$cor[(g-1)*2*q4+(b-1)*q4+1:q4])
                  R <- cov2cor(t(RSr[[b]]) %*% RSr[[b]])
                  indices <- (b-1)*q+1:q         
                  if (Homoc==TRUE) {
                    Sigma[indices,indices] <- R * outer(Stdev[indices],Stdev[indices])
                  } else {
                    Sigmak[indices,indices,g] <- R * outer(Stdev[indices],Stdev[indices])
                  }
                }
              }
            }  
          }
        }
        if (Homoc==TRUE) dimnames(Sigma) <- list(Vnames,Vnames) 
      }  else { 
         Sigma <- ISol$Sigma ; Sigmak <- ISol$Sigmak
       } 
       if (Homoc==TRUE) {
         LnLik <- MClusLikpar(X,n,p,k,tau=tau,muk=muk,Sigma=Sigma,Homoc=TRUE)
       } else {
         LnLik <- MClusLikpar(X,n,p,k,tau=ISol$tau,muk=muk,Sigmak=Sigmak,Homoc=FALSE)
       }  
     }

     Vnames <- names(X)
     Onames <- rownames(X)
     Cnames <- paste("CP",1:k,sep="")
     if (Homoc==TRUE) {
       res <- .Call( "CEMGauss", as.matrix(X), k, Config, TRUE, maxiter, protol, convtol, z,tau,muk,Sigma,NULL,LnLik, FALSE, PACKAGE = "MAINT.Data" )
       rownames(res$muk) <- rownames(res$Sigma) <- colnames(res$Sigma) <- Vnames
       names(res$clusters) <- Onames
       if ( any(!is.finite(res$tau)) || any(!is.finite(res$muk)) || any(!is.finite(res$Sigma)) || CheckSigmaSing(Config,res$Sigma)==TRUE ) 
         res$LnLik <- -Inf  
     }  else {
       res <- .Call( "CEMGauss", as.matrix(X), k, Config, FALSE, maxiter, protol, convtol, z,tau,muk,NULL,Sigmak,LnLik, PACKAGE = "MAINT.Data" )

       rownames(ISol$muk) <- dimnames(ISol$Sigmak)[[1]]  <- dimnames(ISol$Sigmak)[[2]] <- Vnames
       names(res$clusters) <- Onames
       dimnames(res$Sigmak)[[3]] <- Cnames
       if ( any(!is.finite(res$tau)) || any(!is.finite(res$muk)) || any(!is.finite(res$Sigmak)) || CheckSigmakSing(Config,res$Sigmak)==TRUE )
         res$LnLik <- -Inf  
     }

     rownames(res$z) <- Onames
     names(res$tau) <- colnames(res$muk) <- colnames(res$z) <- Cnames

     res$clusters <- Cnames[res$clusters] 
     Lik[rep] <- res$LnLik
     if (res$LnLik > BestLnLik) {
       ISol <- res
       BestLnLik <- res$LnLik 
     }
   }

   return(list(BestSol=ISol,alllnLik=Lik))
} 

