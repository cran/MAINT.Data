setMethod("Idtmclust",
  signature(Idt = "IData"),
  function(Idt, G=1:9, CovCase=1:4, SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","HomandHet"), control=EMControl())   
  {
    call <- match.call()
    SelCrit <- match.arg(SelCrit)
    nrep <- control@nrep
    maxiter <- control@maxiter 
    set.seed(control@seed)
    protol <- control@protol
    convtol <- control@convtol

    X <- data.frame(c(Idt@MidP,Idt@LogR),row.names=rownames(Idt))
    n <- nrow(X)
    p <- ncol(X)
    nCases <- length(CovCase)
    nG <- length(G)
    ntot <- nG*nCases
    GCnames <- paste(rep(paste("G",G,sep=""),each=nCases),paste("C",CovCase,sep=""),sep="")
    Vnames <- names(X)
    Onames <- rownames(X)
    if (Mxt=="HomandHet" || Mxt=="Hom") {
      RepresHom <- vector("list",ntot)
      HomlogLiks <- numeric(ntot)
      HomBICs <- numeric(ntot) 
      HomAICs <- numeric(ntot)
      names(RepresHom) <- names(HomlogLiks) <- names(HomBICs) <- names(HomAICs) <- paste("Hom",GCnames,sep="")
    } else {
      RepresHom <- HomlogLiks <- HomBICs <- HomAICs <- NULL
    }
    if (Mxt=="HomandHet" || Mxt=="Het") {
      RepresHet <- vector("list",ntot)
      HetlogLiks <- numeric(ntot)
      HetBICs <- numeric(ntot) 
      HetAICs <- numeric(ntot)
      names(RepresHet) <- names(HetlogLiks) <- names(HetBICs) <- names(HetAICs) <- paste("Het",GCnames,sep="")
    } else {
      RepresHet <- HetlogLiks <- HetBICs <- HetAICs <- NULL
    }
  
    ind0 <- 0
    for (k in G) {

      if (k==1) {
        for (Cfi in 1:nCases) {
          ind <- ind0+Cfi
          Cf <- CovCase[Cfi]
          Config <- ifelse(Cf==1,1,Cf+1)
          mleEst <- mle(Idt,CovCase=Cf)
          muSig <- coef(mleEst)
          Lik <- mleEst@logLiks[mleEst@BestModel]
          BIC <- mleEst@BICs[mleEst@BestModel]
          AIC <- mleEst@AICs[mleEst@BestModel]
          if (Config==1) { 
            Sigmapar <- p*(p+1)/2 
          }	else if (Config==3) {
            Sigmapar <- 3*p/2	
          }  else if (Config==4) {
            Sigmapar <- p*(p/2+1)/2   
          }	else if (Config==5) {
            Sigmapar <- p
          }  
          if (Mxt=="HomandHet" || Mxt=="Hom") {
            RepresHom[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,Conf=Cf,nG=k,
              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=NULL,classification=NULL,
              parameters=list(pro=NULL,mean=matrix(muSig$mu,ncol=1,dimnames=list(Vnames,NULL)),
                covariance=array(muSig$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL)))
              )
            HomlogLiks[ind] <- Lik
            HomBICs[ind] <- BIC
            HomAICs[ind] <- AIC
          }
          if (Mxt=="HomandHet" || Mxt=="Het") {
            RepresHet[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,Conf=Cf,nG=k,
              logLik=Lik,bic=BIC,aic=AIC,z=NULL,classification=NULL,
              parameters=list(pro=NULL,mean=matrix(muSig$mu,ncol=1,dimnames=list(Vnames,NULL)),
                covariance=array(muSig$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL)))
              )
            HetlogLiks[ind] <- Lik
            HetBICs[ind] <- BIC
            HetAICs[ind] <- AIC
          }
        }  

      } else { 
        ptau <- rep(0.25,k-1)
        sdmuk <- rep(5.,k*p)
        sdSigmaSr <- rep(5.,k*p*(p+1)/2)
        z <- matrix(0.,n,k)
        if (Mxt=="HomandHet" || Mxt=="Hom") {
          clstiallk <- hcEEE(X,minclus=2)
          clsti <- hclass(clstiallk,k)
          for (i in 1:k) z[clsti==i,i] <- 1.
          HomISol <- list(LnLik=-Inf,z=z)
        }  
        if (Mxt=="HomandHet" || Mxt=="Het") {
          clstiallk <- hcVVV(X,minclus=2)
          clsti <- hclass(clstiallk,k)
          z[,] <- 0.
          for (i in 1:k) z[clsti==i,i] <- 1.
          HetISol <- list(LnLik=-Inf,z=z)
        }  

        for (Cfi in 1:nCases) {
          Cf <- CovCase[Cfi]
          Config <- ifelse(Cf==1,1,Cf+1)
          Cf <- CovCase[Cfi]
          Config <- ifelse(Cf==1,1,Cf+1)
          ind <- ind0+Cfi
          if (Mxt=="HomandHet" || Mxt=="Hom") {
            ISol <- EMGaus(X,k,HomISol,maxiter=maxiter,Config=Config,Homoc=TRUE,convtol=convtol,tautol=protol)
            FSol <- RepEMGauss(X,n,p,k,ISol,pertub=list(z=NULL,tau=ptau,muk=sdmuk,SigmaSr=sdSigmaSr),
              nrep=nrep,Config=Config,Homoc=TRUE,maxiter=maxiter,protol=protol,convtol=convtol)
            RepresHom[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,Conf=Cf,nG=k,
              logLik=FSol$BestSol$LnLik,alllnLik=FSol$alllnLik,bic=FSol$BestSol$BIC,aic=FSol$BestSol$AIC,
              parameters=list(
                pro=FSol$BestSol$tau,mean=FSol$BestSol$muk,covariance=array(FSol$BestSol$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL))
              ),
              z=FSol$BestSol$z,classification=FSol$BestSol$clusters
            )
            HomlogLiks[ind] <- FSol$BestSol$LnLik
            HomBICs[ind] <- FSol$BestSol$BIC
            HomAICs[ind] <- FSol$BestSol$AIC
          }
          if (Mxt=="HomandHet" || Mxt=="Het") {
            ISol <- EMGaus(X,k,HetISol,maxiter=maxiter,Config=Config,Homoc=FALSE,convtol=convtol,tautol=protol)
            FSol <- RepEMGauss(X,n,p,k,ISol,pertub=list(z=NULL,tau=ptau,muk=sdmuk,SigmaSr=sdSigmaSr),
              nrep=nrep,Config=Config,Homoc=FALSE,maxiter=maxiter,protol=protol,convtol=convtol)
            RepresHet[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,Conf=Cf,nG=k,
              logLik=FSol$BestSol$LnLik,alllnLik=FSol$alllnLik,bic=FSol$BestSol$BIC,aic=FSol$BestSol$AIC,
              parameters=list(pro=FSol$BestSol$tau,mean=FSol$BestSol$muk,covariance=FSol$BestSol$Sigmak),
              z=FSol$BestSol$z,classification=FSol$BestSol$clusters
            )
            HetlogLiks[ind] <- FSol$BestSol$LnLik
            HetBICs[ind] <- FSol$BestSol$BIC
            HetAICs[ind] <- FSol$BestSol$AIC
          }
        }
      }
      ind0 <- ind0 + nCases
    }  
    if (Mxt=="HomandHet" || Mxt=="Hom") {
      if (SelCrit=="BIC") {
        bestHomMod <- which.min(HomBICs)
        bestHomCrit <- HomBICs[bestHomMod]
      } else if (SelCrit=="AIC") {
        bestHomMod <- which.min(HomAICs)
        bestHomCrit <- HomAICs[bestHomMod]
      }
    } else {
      bestHomMod <- bestHomCrit <- NULL
    }
    if (Mxt=="HomandHet" || Mxt=="Het") {
      if (SelCrit=="BIC") {
        bestHetMod <- which.min(HetBICs)
        bestHetCrit <- HetBICs[bestHetMod]
      } else if (SelCrit=="AIC") {
        bestHetMod <- which.min(HetAICs)
        bestHetCrit <- HetAICs[bestHomMod]
      }
    } else {
      bestHetMod <- bestHetCrit <- NULL
    }

    if ( is.null(bestHetMod) || (!is.null(bestHomMod) && bestHomCrit < bestHetCrit) )
    {
      BestMxt <- "Hom"
    } else {
      BestMxt <- "Het"
    }

     if (BestMxt=="Hom") {
       return (
         new("IdtMclust",call=call,data=Idt,NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,
           BestC=RepresHom[[bestHomMod]]@Conf,BestG=RepresHom[[bestHomMod]]@nG,
           logLiks=c(HomlogLiks,HetlogLiks),BICs=c(HomBICs,HetBICs),AICs=c(HomAICs,HetAICs),
           logLik=RepresHom[[bestHomMod]]@logLik,bic=RepresHom[[bestHomMod]]@bic,aic=RepresHom[[bestHomMod]]@aic,
           parameters=RepresHom[[bestHomMod]]@parameters,z=RepresHom[[bestHomMod]]@z,classification=RepresHom[[bestHomMod]]@classification,
           allres=list(RepresHom=RepresHom,RepresHet=RepresHet) )
         )
     }  else if (BestMxt=="Het") {  
       return (
         new("IdtMclust",call=call,data=Idt,NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,
           BestC=RepresHet[[bestHetMod]]@Conf,BestG=RepresHet[[bestHetMod]]@nG,
           logLiks=c(HomlogLiks,HetlogLiks),BICs=c(HomBICs,HetBICs),AICs=c(HomAICs,HetAICs),
           logLik=RepresHet[[bestHetMod]]@logLik,bic=RepresHet[[bestHetMod]]@bic,aic=RepresHet[[bestHetMod]]@aic,
           parameters=RepresHet[[bestHetMod]]@parameters,z=RepresHet[[bestHetMod]]@z,classification=RepresHet[[bestHetMod]]@classification,
           allres=list(RepresHom=RepresHom,RepresHet=RepresHet) )
         )
     }
  }  
)


EMControl <- function (nrep=100, maxiter=1000, convtol=0.01, protol=1e-6, seed=NULL)  
{
    new("EMControl", nrep=nrep, maxiter=maxiter, convtol=convtol, protol=protol, seed=seed)
}


VectGaussLogLik <- function(x,p,u,SigmaInv,CheckSing=FALSE,cnst=NULL)
{

#  Computes the log-likelihood of the mean, u and covariance matrix, Sigma estimates, 
#  for a multivariate normal distribution, given an  observed vector, x 
#  When CheckSing is set to TRUE it returns a "large penalty" if Sigma is found to be numerically singular
#
#  Parameters:
#  
#  p          -  x dimension
#  x          -  vector of observation data
#  u          -  vector of mean estimates 
#  SigmaInv      -  The inverse of the covariance matrix estimate 
#  CheckSing  -  Boolean flag that should be set to TRUE when it is necessary to check the singularity of Sigma
#  cnst       -  the log-likeklihood constant -1/2 ln(2 Pi |Sigma|), when already known (or null if needs to be computed)

#   value - the log-likelihood of u and Sigma, given x

  if (CheckSing==TRUE)  {
    k2max <- 1E15    # maximum (L2 norm) condition number for Sigma to be considered numerically non-singular 
    PenF <- 1E6
    egvalues <- eigen(SigmaInv,symmetric=TRUE,only.values=TRUE)$values
    k2 <- egvalues[p]/abs(egvalues[1])
    if (!is.finite(k2)) return(-PenF*k2max)
    if (k2 > k2max) return(PenF*(k2max-k2))
  }

  if (is.null(cnst))  {
    lndetSig <- 1./as.double(determinant(SigmaInv,logarithm=TRUE)$modulus)
    cnst <- -0.5 * ( p*log(2*pi) + lndetSig )
  }
  dev <- x-u

  cnst -0.5 * dev %*% SigmaInv %*% dev
}  

MDataGaussLogLik <- function(n,p,X,u,Sigma)
{

#  Computes the log-likelihood of the mean, u and covariance matrix, Sigma estimates, 
#  for a multivariate normal distribution, given an  observed data matrix, X 
#  It returns a "large penalty" if Sigma is found to be numerically singular
#
#  Parameters:
#  
#  p          -  x dimension
#  X          -  vector of observation data
#  u          -  vector of mean estimates 
#  Sigma      -  The covariance matrix estimate 

#   value - a vector with  the log-likelihood contribution of each row of the X data matrix

  k2max <- 1E15      # maximum (L2 norm) condition number for Sigma to be considered numerically non-singular 
  minlndet <- -500   # minimum log-determinant for Sigma to be considered numerically non-singular 
  PenF <- 1E6	   # penalty factor for numerically singular covariance matrices	

  egvalues <- eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values
  k2 <- egvalues[1]/abs(egvalues[p])
  if (!is.finite(k2)) return(rep(-PenF*k2max,n))
  if (k2 > k2max) return(rep(PenF*(k2max-k2),n))

  lndetSig <- as.double(determinant(Sigma,logarithm=TRUE)$modulus)
  if (lndetSig < minlndet) return(rep(PenF*(lndetSig-minlndet),n))
  cnst <- -0.5 * ( p*log(2*pi) + lndetSig )
  SigmaInv <- solve(Sigma) 

  apply(X,1,VectGaussLogLik,p,u,SigmaInv,FALSE,cnst)
}  

EMGaus <- function(X,k,InitSol,Config,Homoc,maxiter,tautol,convtol)
{
   MINZ <- 1E-300
   minzenf <- function(x,defp) if (is.finite(x)) return(max(x,MINZ)) else return(defp)

   n <- nrow(X)
   p <- ncol(X)
   Vnames <- names(X)
   Onames <- rownames(X)
   Wk <- array(0.,dim=c(p,p,k))
   LnLik <- InitSol$LnLik
   Likk <- array(dim=c(n,k))	
   Repmuk <- array(dim=c(n,p,k),dimnames=list(Onames,Vnames,NULL))
   if (is.null(InitSol$z))  {
     z <- array(dim=c(n,k))
     tau <- InitSol$tau
     muk <- InitSol$muk
     if (Homoc==TRUE) {
       Sigma <- InitSol$Sigma
       for (i in 1:k)	Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigma)) 
     } else {
       Sigmak <- InitSol$Sigmak
       for (i in 1:k)	Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigmak[,,i])) 
     }
    Likall <- matrix(rep(apply(Likk,1,sum),k),n,k,byrow=FALSE)
     z <-  matrix(sapply(Likk/Likall,minzenf,defp=1./k),n,k)
   } else {
     z <- InitSol$z
     muk <- matrix(nrow=p,ncol=k,dimnames=list(Vnames,NULL))
     Sigmak <- array(dim=c(p,p,k),dimnames=list(Vnames,Vnames,NULL))
   }
   converg <- FALSE
   iter <- 0
   while (converg == FALSE && iter<maxiter)  {
      nk <- apply(z,2,sum)
      if (any(!is.finite(z)) || any(nk < tautol))  break 	
      else tau <- nk/n	
      for (i in 1:k)  {
        muk[,i] <- apply(matrix(rep(z[,i],p),n,p,byrow=FALSE)*X,2,sum)/nk[i]
        Repmuk[,,i] <- matrix(rep(muk[,i],n),n,p,byrow=TRUE)
	wdev <- as.matrix(matrix(rep(sqrt(z[,i]),p),n,p,byrow=FALSE) * (X - Repmuk[,,i]))
	allind <- 1:p
        midpind <- allind[allind%/%2!=allind/2] 
        lrind <- allind[-midpind]
        if (Config==1) Wk[,,i] <- t(wdev)%*%wdev
	if (Config==4)  {
          Wk[midpind,midpind,i] <- t(wdev[,midpind]) %*% wdev[,midpind] 
          Wk[lrind,lrind,i] <- t(wdev[,lrind]) %*% wdev[,lrind] 
	}
        if (Config==3 || Config==5 ) Wk[,,i] <- diag(apply(wdev^2,2,sum))
        if (Config==3) for (j in 1:(p/2)) {
	  Wk[2*j,2*j-1,i] <- sum(wdev[,2*j-1]*wdev[,2*j])
          Wk[2*j-1,2*j,i] <- Wk[2*j,2*j-1,i]
	}
        if (Homoc==FALSE) {
          Sigmak[,,i] <- Wk[,,i]/nk[i]
          Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigmak[,,i])) 
        }
      }	
      if (Homoc==TRUE) {
        Sigma <- apply(Wk,c(1,2),sum)/n
        dimnames(Sigma) <- list(Vnames,Vnames) 
	for (i in 1:k)	Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigma)) 
      }	
      Likall <- matrix(rep(apply(Likk,1,sum),k),n,k,byrow=FALSE)
      z <-  matrix(sapply(Likk/Likall,minzenf,defp=1./k),n,k)
      LnLik1 <- sum(log(Likall[,1]))
      if (!is.finite(LnLik1)) { LnLik1 <- -Inf ; converg <- TRUE  } 
      if (is.finite(LnLik) && is.finite(LnLik1)) if (abs(LnLik-LnLik1) < convtol) converg <- TRUE 
      if (is.finite(LnLik1)) LnLik <- LnLik1
      iter <- iter+1
    }
    grp <- apply(z,1,which.max)
    names(grp) <- rownames(X)
    if (Config==1) Sigmapar <- p*(p+1)/2	
    else if (Config==3) Sigmapar <- 3*p/2	
    else if (Config==4) Sigmapar <- p*(p/2+1)/2	
    else if (Config==5) Sigmapar <- p
    if (Homoc==TRUE) npar <- p*k + Sigmapar + k-1		
    else npar <- (p+Sigmapar)*k + k-1
    BIC <- -2*LnLik + log(n)*npar	 		
    AIC <- -2*LnLik + 2*npar	 		
    if (Homoc==TRUE) return(list(tau=tau,muk=muk,Sigma=Sigma,z=z,clusters=grp,LnLik=LnLik,npar=npar,BIC=BIC,AIC=AIC))
    else return(list(tau=tau,muk=muk,Sigmak=Sigmak,z=z,clusters=grp,LnLik=LnLik,npar=npar,BIC=BIC,AIC=AIC))
}

MClusLikz <- function(z,X,n,p,k,Config,Homoc,penalty=-1.E99)
{
  Vnames <- names(X)
  Onames <- rownames(X)
  Likk <- array(dim=c(n,k))
  muk <- array(dim=c(p,k),dimnames=list(Vnames,NULL))	
  Repmuk <- array(dim=c(n,p,k),dimnames=list(Onames,Vnames,NULL))
  Wk <- array(dim=c(p,p,k))	
  if (Homoc==FALSE) Sigmak <- array(dim=c(p,p,k))	

  for (i in 1:k)  {
    nk <- apply(z,2,sum)
    tau <- nk/n
    muk[,i] <- apply(matrix(rep(z[,i],p),n,p,byrow=FALSE)*X,2,sum)/nk[i]
    Repmuk[,,i] <- matrix(rep(muk[,i],n),n,p,byrow=TRUE)
    wdev <- as.matrix(matrix(rep(sqrt(z[,i]),p),n,p) * (X - Repmuk[,,i]))
    Wk[,,i] <- t(wdev)%*%wdev
    if (Homoc==FALSE) {
      Sigmak[,,i] <- Wk[,,i]/nk[i]
      Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigmak[,,i])) 
    }
  }	
  if (Homoc==TRUE) {
    Sigma <- apply(Wk,c(1,2),sum)/n
    dimnames(Sigma) <- list(Vnames,Vnames) 
    for (i in 1:k)	Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigma)) 
  }
  lnLik <- log(apply(Likk,1,sum))
  if (any(!is.finite(lnLik))) {
     return(penalty)
  } else {
    return(sum(lnLik))
  }
}

MClusLikpar <- function(X,n,p,k,tau,muk,Sigma=NULL,Sigmak=NULL,Homoc,penalty=-1.E99)
{
  Likk <- array(dim=c(n,k))

  for (i in 1:k)  {
   	 if (Homoc==TRUE) Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigma)) 
         else Likk[,i] <- tau[i] * exp(MDataGaussLogLik(n,p,X,muk[,i],Sigmak[,,i])) 
   }
   lnLik <- log(apply(Likk,1,sum))
   if (any(!is.finite(lnLik))) return(penalty)
   else return(sum(lnLik))
}

RepEMGauss <- function(X,n,p,k,ISol,pertub,nrep,Config,Homoc,maxiter,convtol,protol)
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
      if (!is.null(pertub$SigmaSr)) {
        q0 <- p*(p-1)/2
       	if (Homoc==TRUE) {
       	  Sk <- 1
       	} else {
       	  Sk <- k
       	}  
        for (i in 1:Sk)  {
          if (Homoc==TRUE) {
            SigmaSr <- mychol(ISol$Sigma,p)
          } else {
            SigmaSr <- mychol(ISol$Sigmak[,,i],p)
          } 
          SigmaSr[col(SigmaSr)>row(SigmaSr)] <- 
            rnorm(q0,SigmaSr[col(SigmaSr)>row(SigmaSr)],pertub$SigmaSr[((i-1)*q0+1):(i*q0)])
          SigmaSr[col(SigmaSr)==row(SigmaSr)] <- 
            abs(rnorm(p,SigmaSr[col(SigmaSr)==row(SigmaSr)],pertub$SigmaSr[(Sk*q0+(i-1)*p+1):(Sk*q0+i*p)]))
          if (Homoc==TRUE) {
            Sigma <- t(SigmaSr) %*% SigmaSr
            dimnames(Sigma) <- list(Vnames,Vnames) 
          } else {
            Sigmak[,,i] <- t(SigmaSr) %*% SigmaSr
          } 
         }
       }  else { 
         Sigma <- ISol$Sigma ; Sigmak <- ISol$Sigmak
       } 
       if (Homoc==TRUE) {
         LnLik <- MClusLikpar(X,n,p,k,tau=tau,muk=muk,Sigma=Sigma,Homoc=TRUE)
       } else {
         LnLik <- MClusLikpar(X,n,p,k,tau=ISol$tau,muk=muk,Sigmak=Sigmak,Homoc=FALSE)
       }  
     }

     if (Homoc==TRUE) {
       res <- EMGaus(X,k,list(z=z,tau=tau,muk=muk,Sigma=Sigma,LnLik=LnLik),Config=Config,Homoc=TRUE,maxiter=maxiter,convtol=convtol,tautol=protol)
     }  else {
       res <- EMGaus(X,k,list(z=z,tau=tau,muk=muk,Sigmak=Sigmak,LnLik=LnLik),Config=Config,Homoc=FALSE,maxiter=maxiter,convtol=convtol,tautol=protol)
     }  
     Lik[rep] <- res$LnLik
     if (res$LnLik > BestLnLik) {
       ISol <- res
       BestLnLik <- res$LnLik 
     }
   }

   return(list(BestSol=ISol,alllnLik=Lik))
} 


