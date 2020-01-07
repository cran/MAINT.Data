setMethod("Idtmclust",
  signature(Idt = "IData"),
  function(Idt, G=1:9, CovCase=1:4, SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","HomandHet"), control=EMControl())   
  {
    pertubzVct <- function(z,maxprt) {
      maxzind <- which.max(z)
      change <- runif(1,0.,min(z[maxzind],maxprt))
      z[maxzind] <- z[maxzind] - change
      z[-maxzind] <- z[-maxzind] + change/(length(z)-1)
      z
    } 
    pertubzMat <- function(OrigzMat,maxprt=0.1) t(apply(OrigzMat,1,pertubzVct,maxprt=maxprt))

    call <- match.call()
    SelCrit <- match.arg(SelCrit)
    Mxt <- match.arg(Mxt)
    nrep <- control@nrep
    maxiter <- control@maxiter 
    set.seed(control@seed)
    protol <- control@protol
    convtol <- control@convtol

    if (class(G)!="integer") G <- as.integer(G)
    if (class(CovCase)!="integer") CovCase <- as.integer(CovCase)

    X <- data.frame(c(Idt@MidP,Idt@LogR),row.names=rownames(Idt))
    n <- nrow(X)
    p <- ncol(X)
    q <- p/2
    if (q==1) CovCase <- q1CovCase(CovCase) 
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
            clusters <- rep("CP1",Idt@NObs)
            names(clusters) <- Idt@ObsNames
            RepresHom[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,Conf=Cf,nG=k,
#              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=NULL,classification=rep("CP1",Idt@NObs),
              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=NULL,classification=clusters,
              parameters=list(pro=NULL,mean=matrix(muSig$mu,ncol=1,dimnames=list(Vnames,"CP1")),
                covariance=array(muSig$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL)))
              )
            HomlogLiks[ind] <- Lik
            HomBICs[ind] <- BIC
            HomAICs[ind] <- AIC
          }
          if (Mxt=="HomandHet" || Mxt=="Het") {
            clusters <- rep("CP1",Idt@NObs)
            names(clusters) <- Idt@ObsNames
            RepresHet[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,Conf=Cf,nG=k,
#              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=NULL,classification=rep("CP1",Idt@NObs),
              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=NULL,classification=clusters,
              parameters=list(pro=NULL,mean=matrix(muSig$mu,ncol=1,dimnames=list(Vnames,"CP1")),
                covariance=array(muSig$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,"CP1")))
              )
            HetlogLiks[ind] <- Lik
            HetBICs[ind] <- BIC
            HetAICs[ind] <- AIC
          }
        }  

      } else { 
        pertubfct <- 2.                               # pertubation factors
        ptaufct <- 0.2   
        pR <- 0.1
        ptau <- pertubfct*rep(ptaufct/k,k-1)
        
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

        Vnames <- names(X)
        Onames <- rownames(X)
        Cnames <- paste("CP",1:k,sep="")

        for (Cfi in 1:nCases) {
          Cf <- CovCase[Cfi]
          Config <- ifelse(Cf==1,1,Cf+1)
          Cf <- CovCase[Cfi]
          Config <- ifelse(Cf==1,1,Cf+1)
          ind <- ind0+Cfi

          ISol <- NULL 
          maxzpert <- 0.1
          maxinittrials <- 100
          inittrial <- 0
          if (Mxt=="HomandHet" || Mxt=="Hom") {
            while ( is.null(ISol)  || any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigma)) )  
            {
              inittrial <- inittrial + 1
              ISol <- .Call( "CEMGauss", as.matrix(X), k, Config, TRUE, maxiter, protol, convtol,
#                HomISol$z, NULL, NULL, NULL, NULL, HomISol$LnLik, 
                 HomISol$z, NULL, NULL, NULL, NULL, HomISol$LnLik, FALSE, 
                 PACKAGE = "MAINT.Data" )
              rownames(ISol$muk) <- rownames(ISol$Sigma) <- colnames(ISol$Sigma) <- Vnames
              names(ISol$clusters) <- Onames
              if ( any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigma)) ) {
                if (inittrial < maxinittrials) HomISol$z <- pertubzMat(HomISol$z)
                else stop(paste("Idtmclust failed to find a valid Homoscedastic ",k,"-group solution for configuration C",CovCase," after ",inittrial,"trials.\n",sep=""))
              }
            } 

            rownames(ISol$z) <- Onames
            names(ISol$tau) <- colnames(ISol$muk) <- colnames(ISol$z) <- Cnames
            ISol$clusters <- Cnames[ISol$clusters] 
            names(ISol$clusters) <- Onames
            stdv <- sqrt(diag(ISol$Sigma))
            sdmuk <- pertubfct*rep(stdv,k)
            sdStdev <- pertubfct*stdv/4
            if (Config==1 || Config==2) sdR <- pertubfct*rep(pR,p*(p-1)/2)
            else if (Config==3) sdR <- pertubfct*rep(pR,p)
            else if (Config==4) sdR <- pertubfct*rep(pR,q*(q-1))
            else if (Config==5) sdR <- NULL
            
            FSol <- RepEMGauss(X,n,p,k,ISol,pertub=list(z=NULL,tau=ptau,muk=sdmuk,Stdev=sdStdev,cor=sdR),
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
            while ( is.null(ISol)  || any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigmak)) )  
            {
              inittrial <- inittrial + 1
                ISol <- .Call( "CEMGauss", as.matrix(X), k, Config, FALSE, maxiter, protol, convtol,
#                  HetISol$z, NULL, NULL, NULL, NULL, HetISol$LnLik, 
                  HetISol$z, NULL, NULL, NULL, NULL, HetISol$LnLik, FALSE, 
                  PACKAGE = "MAINT.Data" )
                rownames(ISol$muk) <- dimnames(ISol$Sigmak)[[1]]  <- dimnames(ISol$Sigmak)[[2]] <- Vnames
                names(ISol$clusters) <- Onames
              if ( any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigmak)) ) {
                if (inittrial < maxinittrials) HetISol$z <- pertubzMat(HetISol$z)
                else stop(paste("Idtmclust failed to find a valid Heteroscedastic ",k,"-group solution for configuration C",CovCase," after ",inittrial,"trials.\n",sep=""))
              }
            }

            dimnames(ISol$Sigmak)[[3]] <- Cnames
            rownames(ISol$z) <- Onames
            names(ISol$tau) <- colnames(ISol$muk) <- colnames(ISol$z) <- Cnames
            ISol$clusters <- Cnames[ISol$clusters] 
            names(ISol$clusters) <- Onames
            q0 <- p*(p-1)/2 
            sdmuk <- numeric(p*k)
            sdStdev <- numeric(p*k)
            for (g in 1:k) {
              stdv <- sqrt(diag(ISol$Sigmak[,,g]))
              sdmuk[(g-1)*p+1:p] <- pertubfct*stdv
              sdStdev[(g-1)*p+1:p] <- pertubfct*stdv/4
            }
            if (Config==1 || Config==2) sdR <- pertubfct*rep(pR,k*p*(p-1)/2)
            else if (Config==3) sdR <- pertubfct*rep(pR,k*p)
            else if (Config==4) sdR <- pertubfct*rep(pR,k*q*(q-1))
            else if (Config==5) sdR <- NULL
            
            FSol <- RepEMGauss(X,n,p,k,ISol,pertub=list(z=NULL,tau=ptau,muk=sdmuk,Stdev=sdStdev,cor=sdR),
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
           logLiks=c(HetlogLiks,HetlogLiks),BICs=c(HetBICs,HetBICs),AICs=c(HetAICs,HetAICs),
           logLik=RepresHet[[bestHetMod]]@logLik,bic=RepresHet[[bestHetMod]]@bic,aic=RepresHet[[bestHetMod]]@aic,
           parameters=RepresHet[[bestHetMod]]@parameters,z=RepresHet[[bestHetMod]]@z,classification=RepresHet[[bestHetMod]]@classification,
           allres=list(RepresHet=RepresHet,RepresHet=RepresHet) )
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
#  X          -  matrix of observation data
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
   LikExp <- matrix(nrow=n,ncol=k)

   Repmuk <- array(dim=c(n,p,k),dimnames=list(Onames,Vnames,NULL))
   if (is.null(InitSol$z))  {
     z <- array(dim=c(n,k))
     tau <- InitSol$tau
     muk <- InitSol$muk
     if (Homoc==TRUE) {
       Sigma <- InitSol$Sigma
        for (i in 1:k)	LikExp[,i] <- MDataGaussLogLik(n,p,X,muk[,i],Sigma) 
     } else {
       Sigmak <- InitSol$Sigmak
       for (i in 1:k)	LikExp[,i] <- MDataGaussLogLik(n,p,X,muk[,i],Sigmak[,,i]) 
     }
     nrmfct <- apply(LikExp,1,max)
     Likk <- scale(exp(sweep(LikExp,1,nrmfct)),center=FALSE,scale=1/tau)
     attr(Likk,"scaled:scale") <- NULL
     Likall <- rowSums(Likk)
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
        q <- p/2
        midpind <- allind[1:q] 
        lrind <- allind[-midpind]
        if (Config==1) Wk[,,i] <- t(wdev)%*%wdev
	      if (Config==4)  {
          Wk[midpind,midpind,i] <- t(wdev[,midpind]) %*% wdev[,midpind] 
          Wk[lrind,lrind,i] <- t(wdev[,lrind]) %*% wdev[,lrind] 
	      }
        if (Config==3 || Config==5 ) Wk[,,i] <- diag(apply(wdev^2,2,sum))
        if (Config==3) for (j in 1:(p/2)) {
	        Wk[q+j,j,i] <- Wk[j,q+j,i] <- sum(wdev[,j]*wdev[,q+j])
	      }
        if (Homoc==FALSE) {
          Sigmak[,,i] <- Wk[,,i]/nk[i]
          LikExp[,i] <- MDataGaussLogLik(n,p,X,muk[,i],Sigmak[,,i]) 
        }
      }	

      if (Homoc==TRUE) {
        Sigma <- apply(Wk,c(1,2),sum)/n
        dimnames(Sigma) <- list(Vnames,Vnames) 
        for (i in 1:k)	LikExp[,i] <- MDataGaussLogLik(n,p,X,muk[,i],Sigma) 
      }	
      nrmfct <- apply(LikExp,1,max)
      Likk <- scale(exp(sweep(LikExp,1,nrmfct)),center=FALSE,scale=1/tau)
      attr(Likk,"scaled:scale") <- NULL
      Likall <- rowSums(Likk)
      z <-  matrix(sapply(Likk/Likall,minzenf,defp=1./k),n,k)
      LnLik1 <- sum(nrmfct+log(Likall))

      if (!is.finite(LnLik1)) { LnLik <- -Inf ; converg <- TRUE  } 
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



