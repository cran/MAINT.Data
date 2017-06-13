SNCnf1MaxLik <- function(Data,grouping=NULL,initpar=NULL,EPS=1E-6,OptCntrl=list())
{
#  Note  -  The Data argument should be a matrix containing the mid-points in the first columns and the log-ranges in the next columns 

  n <- nrow(Data)		# number of observations
  p <- ncol(Data)		# total number of variables (mid-points + log-ranges)
  q <- p/2			# number of interval variables

  sd0_default <- 0.1
  if (!is.null(OptCntrl$sd0))  {
    sd0 <- OptCntrl$sd0
  }  else  {
    sd0 <- sd0_default
  } 
  OptCntrl$sd0 <- NULL 
  maxiter <- 500
  maxeval <- 1000
  reltol <- sqrt(.Machine$double.eps)
  mymsncontrol <- list(iter.max=maxiter,eval.max=maxeval,rel.tol=reltol)
  OptCntrl$iter.max <- maxiter
  OptCntrl$eval.max <- maxeval
  OptCntrl$rel.tol <- reltol

  if (is.null(grouping))  {
    k <- 1  
    parsd <- rep(sd0,2*p)      	# standard deviation hyper-parameters - used to generate random starting points
    if (is.null(initpar))  { initpar <- mymsn.mle(NULL,y=Data,control=mymsncontrol)$par }
    res <- RepLOptim(initpar,parsd,fr=NULL,method=mymsn.mle,y=Data,control=OptCntrl)
    ksi <- res$par[1:p]
    beta2k <- NULL
    dev <- scale(Data,center=ksi,scale=FALSE)
  }  else {
    lev <- levels(grouping)
    k <- length(lev)
    parsd <- rep(sd0,(k+1)*p)      	# standard deviation hyper-parameters - used to generate random starting points
    X <- model.matrix(Data ~ grouping)
    if (is.null(initpar))  { initpar <- mymsn.mle(NULL,y=Data,x=X,control=mymsncontrol)$par }
    res <- RepLOptim(initpar,parsd,fr=NULL,method=mymsn.mle,y=Data,x=X,control=OptCntrl)
    beta <- matrix(res$par[1:(k*p)],nrow=k,ncol=p)
    beta2k <- beta[-1,]
    ksi <- matrix(nrow=k,ncol=p)
    ksi[1,] <- beta[1,]
    ksi[2:k,] <- scale(beta2k,center=-ksi[1,],scale=FALSE)
    dev <- matrix(nrow=n,ncol=p)
    for (g in 1:k)  {
     gind <- which(grouping==lev[g])
     dev[gind,] <- scale(Data[gind,],center=ksi[g,],scale=FALSE)
    }
  }
  Omega <- t(dev) %*% dev/n
  eta  <- res$par[(k*p+1):((k+1)*p)]
  alpha <- eta*sqrt(diag(Omega))
  Omegabar <- cov2cor(Omega)
  tmp <- drop(Omegabar%*%alpha)
  c2 <- 1. + alpha%*%tmp
#  delta <- tmp/sqrt(c2)
  delta <- tmp/rep(sqrt(c2),length(tmp))
  if (is.null(grouping))
  {
    CP <- cnvDPtoCP(p,ksi,Omega,alpha)
    mu <- CP$mu
  }  else {
    mu <- matrix(nrow=k,ncol=p)
    CP <- cnvDPtoCP(p,ksi[1,],Omega,alpha)
    mu[1,] <- CP$mu
    mu[2:k,] <- scale(beta2k,center=-mu[1,],scale=FALSE)
  }
  list(lnLik=res$val/(-2),ksi=ksi,beta2k=beta2k,Omega=Omega,Omega.cor=Omegabar,alpha=alpha,
    delta=delta,mu=mu,Sigma=CP$Sigma,gamma1=CP$gamma1,admissible=TRUE,c2=c2,optres=res)
}

Getvarofsmplmeansbygrp <- function(X,grouping)
{
  varbygrp <- apply(X,2,function(x) tapply(x,grouping,var))
  nk <- as.numeric(table(grouping))
  sweep(varbygrp,1,nk,"/")
}

FreePar <- function(q,j1,j2,Conf)
{
  if (Conf==1 || j1==j2)  {
    return(TRUE)
  }  
  if (Conf==5) { return(FALSE) }
  if ( ((j1<=q) && (j2<=q)) || ((j1>q) && (j2>q)) )
  {
    if (Conf==2 || Conf==4)  {
      return(TRUE)
    }  else {
      return(FALSE)
    }
  }  else {
    if (j2==j1+q)
    {
      if (Conf==2 || Conf==3)  {
        return(TRUE)
      }  else  {
        return(FALSE)
      }
    }  else {
      return(FALSE)
    }
  }
}  		

SKnpar <- function(Conf,p,q,Ngrps=1,Mxt=c("LocMod","GenMod"))
{
  Mxt <- match.arg(Mxt)
  ksipar <- Ngrps*p 
  if (Mxt=="LocMod")
  {
    alphapar <- p 
    if (Conf==1)  { Omgpar <- p*(p+1)/2 }
    if (Conf==2)  { Omgpar <- q*(q+1) + p }
    if (Conf==3)  { Omgpar <- p + q }
    if (Conf==4)  { Omgpar <- q*(q+1) }
    if (Conf==5)  { Omgpar <- p }
  }  else if (Mxt=="GenMod")  {
    alphapar <- Ngrps*p 
    if (Conf==1)  { Omgpar <- Ngrps*p*(p+1)/2 }
    if (Conf==2)  { Omgpar <- Ngrps*(q*(q+1)+p) }
    if (Conf==3)  { Omgpar <- Ngrps*(p+q) }
    if (Conf==4)  { Omgpar <- Ngrps*q*(q+1) }
    if (Conf==5)  { Omgpar <- Ngrps*p }
  }
  ksipar+alphapar+Omgpar  # return(ksipar+alphapar+Omgpar)
}

SNVCovscaling <- function(Conf,p,stdv,k=1)  # Creates a scaling matrix in order to convert standardized
                                            # vcov matrix, to a vcov matrix in tbe original scale, for the Skew-Normal model    
#  Arguments:

#  Conf    -  Configuration (1,2,3,4 or 5) of variance-covariance matrix 
#  p       -  Number of variables
#  stdv    -  Vector of standard deviations
#  k       -  Number of groups

#  Value     - A matrix with the scaling constants   
{
  if ( p!=length(stdv) ) {
    stop("VCovscaling arguments with wrong dimensions\n")
  }
  nlocpar <- k*p
  npar <- nlocpar + p*(p+1)/2 + p
  sclvct <- numeric(npar)
  sclvct[1:nlocpar] <- rep(stdv,k)
  gamma1ind <- (npar-p+1):npar
  sclvct[gamma1ind] <- rep(1,p)
  ind <- nlocpar
  for(i in 1:p) for(j in i:p)  {
    ind <- ind + 1  
    sclvct[ind] <- stdv[i]*stdv[j]
  }
  if (Conf!=1) { sclvct <- sclvct[c(1:nlocpar,nlocpar+vcovCind(Conf,p),gamma1ind)] } 
  outer(sclvct,sclvct)  
}

IdtSNmle <- function(Idt, grouping=NULL, Type=c("SingDst","HomMxt"), CVtol=1.0e-5, bordertol=1e-2,
  OptCntrl=list(), onerror=c("stop","warning","silentNull"), CovCaseArg, Config, SelCrit)
{

##########################################################################
### Auxiliary functions
##########################################################################

  ImprvSNmdl <- function(RestModel,RestConf,FullConf)
  {
    if (FullConf==1)
    {
      if (Type=="SingDst")  {
        DP <- cnvCPtoDP(p,RestModel$mu,RestModel$Sigma,RestModel$gamma1) 
      }  else if (Type=="HomMxt")  {  
        DP <- cnvCPtoDP(p,RestModel$mu[1,]-RestModel$mu0[1,],RestModel$Sigma,RestModel$gamma1)
      }
      newpar <- c(DP$ksi,DP$alpha/DP$omega)
      FullModel <- SNCnf1MaxLik(Xscld,initpar=newpar,OptCntrl=OptCntrl)
    }  else  {
      SigmaSrpar <- GetCovPar(RestModel$Sigma,FullConf) 
      if (Type=="SingDst")  {
        newpar <- c(RestModel$mu,SigmaSrpar,RestModel$gamma1)
      }  else if (Type=="HomMxt")  {
        newpar <- c(RestModel$mu[1,]-RestModel$mu0[1,],SigmaSrpar,RestModel$gamma1)
      }
      FullModel <- SNCMaxLik(Xscld,Config=FullConf,initpar=newpar,OptCntrl=OptCntrl)
    }
    if (Type=="HomMxt" && !is.null(FullModel$ksi))
    { 
      FullModel$mures <- FullModel$mu
      FullModel$ksires <- FullModel$ksi
      FullModel$mu <- scale(RestModel$mu0,center=-FullModel$mu,scale=FALSE)
      FullModel$ksi <- scale(RestModel$mu0,center=-FullModel$ksi,scale=FALSE)
      attr(FullModel$mu,"scaled:center") <- attr(FullModel$ksi,"scaled:center") <- NULL
    }

    FullModel
  }

  Getvcov <- function(Res)
  {
    if (Type=="SingDst")
    {
      SngDparnam <- paste("mu_",Xnames,sep="")
      for (i in 1:p)
        SngDparnam <- c(SngDparnam,paste("Sigma_",Xnames[i],"_",Xnames[i:p],sep=""))
      SngDparnam <- c(SngDparnam,paste("gamma1_",Xnames,sep=""))
      npar <- SKnpar(Conf,p,p/2)
      InFData <- mysn.infoMv( dp=list(xi=Res$ksi,Omega=Res$Omega,alpha=Res$alpha),
        y=Xscld, w=n1scvct )
      if (InFData$status!="Regular") {
        return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL,status=InFData$status) )
      }
      if (Conf==1)  {
        mleCPvcov <- InFData$asyvar.cp
        rownames(mleCPvcov) <- colnames(mleCPvcov) <- SngDparnam
      }  else  {
        nparC1 <- 2*p + p*(p+1)/2  
        parind <- c(1:p,p+SigCind(Conf,p/2),(nparC1-p+1):nparC1)
        mleCPvcov <- pdwt.solve(InFData$info.cp[parind,parind])
        if ( !is.null(mleCPvcov) ) { rownames(mleCPvcov) <- colnames(mleCPvcov) <- SngDparnam[parind] }
      }
    }  else if (Type=="HomMxt") {
      HoMxtparnam <- paste("mu_",Xnames,"_",grplvls[1],sep="")
      for (g in 2:k) HoMxtparnam <- c(HoMxtparnam,paste("mu_",Xnames,"_",grplvls[g],sep=""))
      for (i in 1:p) HoMxtparnam <- c(HoMxtparnam,paste("Sigma_",Xnames[i],"_",Xnames[i:p],sep=""))
      HoMxtparnam <- c(HoMxtparnam,paste("gamma1_",Xnames,sep=""))
      npar <- SKnpar(Conf,p,p/2,Ngrps=k)
      InFData <- mysn.infoMv( dp=list(beta=rbind(Res$ksi[1,],Res$beta2k),Omega=Res$Omega,alpha=Res$alpha),
        y=Xscld, x=grpModMat, w=n1scvct )
      if (InFData$status!="Regular") {
        return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL,status=InFData$status) )
      }
      if (Conf==1)  {
        betavcov <- InFData$asyvar.cp
        nSpar <- p*(p+1)/2
      }  else  {
        nparC1 <- (k+1)*p + p*(p+1)/2  
        parind <- c(1:(k*p),(k*p)+SigCind(Conf,p/2),(nparC1-p+1):nparC1)
        nSpar <- length(parind) - (k+1)*p
        betavcov <- pdwt.solve(InFData$info.cp[parind,parind])
      }
      if ( !is.null(betavcov) )
      {
        mleCPvcov <- matrix(nrow=npar,ncol=npar)		
        muind <- 1:(p*k)
        nmugind <- NULL
        for (g in 1:k) nmugind <- c(nmugind,(0:(p-1))*k+g)
        Sind <- p*k + 1:nSpar
        gamma1ind <- p*k + nSpar + 1:p
        Sgamma1ind <- c(Sind,gamma1ind)
        M <- kronecker(diag(p),Gmat)
        mleCPvcov[muind,muind] <- (M %*% betavcov[muind,muind] %*% t(M))[nmugind,nmugind]
        mleCPvcov[muind,Sgamma1ind] <- (M %*% betavcov[muind,Sgamma1ind])[nmugind,]
        mleCPvcov[Sgamma1ind,muind] <- t(mleCPvcov[muind,Sgamma1ind])
        mleCPvcov[Sgamma1ind,Sgamma1ind] <- betavcov[Sgamma1ind,Sgamma1ind]
        if (Conf==1)  {
          rownames(mleCPvcov) <- colnames(mleCPvcov) <- HoMxtparnam
        } else  {
          rownames(mleCPvcov) <- colnames(mleCPvcov) <- HoMxtparnam[parind]
        }
      } else { 
        mleCPvcov <- NULL
      }
    }
    if (is.null(mleCPvcov))  {
     return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL) )
    }  
    mleCPvcov <- mleCPvcov * SNVCovscaling(Conf,p,Xsd,k=k) # scale back vcov matrix   
    CPStderr <- sqrt(diag(mleCPvcov))
    if (Type=="SingDst") {
      muEse <- CPStderr[1:p]
    }  else if (Type=="HomMxt") {
      muEse <- matrix(CPStderr[1:(k*p)],nrow=k,ncol=p,byrow=TRUE)
    }
    gammaind <- (npar-p+1):npar
    gamma1Ese <- CPStderr[gammaind]
    SigmaEse <- matrix(nrow=p,ncol=p) 
    cnt <- k*p
    for (j1 in 1:p) for (j2 in j1:p)
    {   
      if (FreePar(q,j1,j2,Conf))
      {
        cnt <- cnt+1
        SigmaEse[j1,j2] <- SigmaEse[j2,j1] <- CPStderr[cnt]
      }
    }
    names(gamma1Ese) <- rownames(SigmaEse) <- colnames(SigmaEse) <- Xnames
    if (Type=="SingDst")  {
      names(muEse) <- Xnames
    }  else if (Type=="HomMxt") {
      rownames(muEse) <- grplvls
      colnames(muEse) <- Xnames
    }

    list(mleCPvcov=mleCPvcov,muEse=muEse,SigmaEse=SigmaEse,gamma1Ese=gamma1Ese,status="Regular")
  }

##########################################################################
### End of auxiliary functions
### Start of main IdtSNmle function
##########################################################################

  Type <- match.arg(Type)
  onerror <- match.arg(onerror)

  q <- Idt@NIVar
  p <- 2*q
  nvcovpar <- p*(p+1)/2
  n <- Idt@NObs
  n1scvct <- rep(1,n)
  maxsk <- 0.99527

  if (Type=="SingDst")  { 
    k <- 1
  }  else if (Type=="HomMxt")  {
    if (is.null(grouping))  { error(onerror,"Argument grouping is missing from MANOVA method\n") }
    if (!is.factor(grouping))  { error(onerror,"'grouping' is not a factor") }
    nk <- as.numeric(table(grouping))
    if (n != sum(nk))  {
      error(onerror,"Number of observations in IData object and grouping factor do not agree with each other\n")
    }
    k <- length(nk) 
    if (k==1)  {
      error(onerror,"The data belongs to one single group. A partition into at least two different groups is required\n")
    }
    grplvls <- levels(grouping)
    Gmat <- model.matrix(~ grplvls)
    grpModMat <- model.matrix(~ grouping)
  }
  Config <- sort(Config)  

  X <- cbind(Idt@MidP,Idt@LogR)
  Xnames <- names(X)
  if (CovCaseArg)  {
    nCovCases <- 4 
    CovCaseMap <- c(1,NA,2,3,4)
    slotnames <- paste("SNModCovC",1:4,sep="")
  } else {
    nCovCases <- 5 
    CovCaseMap <- 1:5
    slotnames <- paste("SNC",1:5,sep="")
  }
  logLiks <- AICs <- BICs <- rep(NA_real_,nCovCases)
  CovConfCases <- vector("list",nCovCases)
  names(logLiks) <- names(AICs) <- names(BICs) <- names(CovConfCases) <- slotnames
  for (model in Config)
  {
    CovConfCases[[CovCaseMap[model]]] <- 
      list( muE=NULL,muEse=NULL,gamma1E=NULL,gamma1Ese=NULL,SigmaE=NULL,SigmaEse=NULL,status=NULL,
        ksiE=NULL,alphaE=NULL,OmegaE=NULL,mleCPvcov=NULL,logLik=NULL,AIC=NULL,BIC=NULL,optres=NULL )
  }
  Xscld <- scale(X)		        # standardized data  
  Xmean <- attr(Xscld,"scaled:center")  
  Xsd <- attr(Xscld,"scaled:scale")  
  lglikdif <- -n*sum(log(Xsd))		# difference between the log-likelihoods for the original and standardized data  
  StdDtRes <- vector("list",nCovCases)
  noestmsg <- "No model was estimated for Covariance configuration case"
  for (i in Config) StdDtRes[[CovCaseMap[i]]] <- paste(noestmsg,CovCaseMap[i])
  StdDtRes4 <- NULL 
  for (Conf in Config)
  {
    if (Conf==1)  {
      prevMod <- NULL
    }  else if (Conf == 2)  {
      if (is.character(StdDtRes[[1]])) { 
        prevMod <- NULL
      }  else  {
        prevMod <- StdDtRes[[1]]
      }
    }  else if (Conf==3 || Conf==4)
    { 
      if (CovCaseArg || is.character(StdDtRes[[2]]))
      { 
        if (is.character(StdDtRes[[1]]))  {
          prevMod <- NULL
        }  else  {
          prevMod <- StdDtRes[[1]]
        }
      }  else  {
        prevMod <- StdDtRes[[2]]
      }
    }  else if (Conf==5)  {
      if (is.character(StdDtRes[[CovCaseMap[3]]]))
      { 
        if (is.character(StdDtRes[[CovCaseMap[4]]]))
        { 
          if (CovCaseArg || is.character(StdDtRes[[2]]))
          { 
            if (is.character(StdDtRes[[1]]))  {
              prevMod <- NULL
            }  else {
              prevMod <- StdDtRes[[1]]
            }
          }  else { 
            prevMod <- StdDtRes[[2]]
          }
        } else {
          prevMod <- StdDtRes[[CovCaseMap[4]]]
        }
      }  else {
        prevMod <- StdDtRes[[CovCaseMap[3]]]
      }
    }
    if (Type=="SingDst")
    {
      if (Conf == 1)  {
        StdDtRes[[1]] <- SNCnf1MaxLik(Xscld,OptCntrl=OptCntrl)
      }  else if (Conf == 2) {
        StdDtRes4 <- SNCMaxLik(Xscld,Config=4,OptCntrl=OptCntrl,prevMod=prevMod)
        StdDtRes[[2]] <- ImprvSNmdl(RestModel=StdDtRes4,RestConf=4,FullConf=2)
      }  else if (Conf == 4 && !is.null(StdDtRes4)) {
        StdDtRes[[CovCaseMap[4]]] <- StdDtRes4
      }  else {
        StdDtRes[[CovCaseMap[Conf]]] <- SNCMaxLik(Xscld,Config=Conf,OptCntrl=OptCntrl,prevMod=prevMod)
      } 
    }  else if (Type=="HomMxt")
    {
      if (Conf == 1)  {
        StdDtRes[[1]] <- SNCnf1MaxLik(Xscld,grouping,OptCntrl=OptCntrl)
      }  else if (Conf == 2) {
        StdDtRes4 <- SNCMaxLik(Xscld,Config=4,grouping=grouping,OptCntrl=OptCntrl,prevMod=prevMod)
        StdDtRes[[2]] <- ImprvSNmdl(RestModel=StdDtRes4,RestConf=4,FullConf=2)
      }  else if (Conf == 4 && !is.null(StdDtRes4))  {
        StdDtRes[[CovCaseMap[4]]] <- StdDtRes4 
      }  else  {
        StdDtRes[[CovCaseMap[Conf]]] <- SNCMaxLik(Xscld,Config=Conf,grouping=grouping,OptCntrl=OptCntrl,prevMod=prevMod) 
      }
    }
  }	

  if ( is.element(4,Config) && is.element(5,Config) )
  { 
    if ( StdDtRes[[CovCaseMap[4]]]$lnLik < StdDtRes[[CovCaseMap[5]]]$lnLik )  {
      StdDtRes[[CovCaseMap[4]]] <- ImprvSNmdl(RestModel=StdDtRes[[CovCaseMap[5]]],RestConf=5,FullConf=4)
    }
  }
  if ( is.element(3,Config) && is.element(5,Config) )
  { 
    if ( StdDtRes[[CovCaseMap[3]]]$lnLik < StdDtRes[[CovCaseMap[5]]]$lnLik )  {
      StdDtRes[[CovCaseMap[3]]] <- ImprvSNmdl(RestModel=StdDtRes[[CovCaseMap[5]]],RestConf=5,FullConf=3)
    }
  }
  if ( is.element(2,Config) && is.element(4,Config) )
  { 
    if ( StdDtRes[[2]]$lnLik < StdDtRes[[CovCaseMap[4]]]$lnLik )  {
      StdDtRes[[2]] <- ImprvSNmdl(RestModel=StdDtRes[[CovCaseMap[4]]],RestConf=4,FullConf=2)
    }
  }
  if ( is.element(2,Config) && is.element(3,Config) )
  { 
    if ( StdDtRes[[2]]$lnLik < StdDtRes[[CovCaseMap[3]]]$lnLik )  {
      StdDtRes[[2]] <- ImprvSNmdl(RestModel=StdDtRes[[CovCaseMap[3]]],RestConf=3,FullConf=2)
    }
  }
  if ( is.element(1,Config) && is.element(2,Config) )  { 
    if ( StdDtRes[[1]]$lnLik < StdDtRes[[2]]$lnLik )  {
      StdDtRes[[1]] <- ImprvSNmdl(RestModel=StdDtRes[[2]],RestConf=2,FullConf=1)
    }
  } 

  for (Conf in Config)
  {
    if (!is.null(StdDtRes[[CovCaseMap[Conf]]]$optres$par))
    {
      CvCase <- CovCaseMap[Conf]	
      if (Type=="SingDst")
      {
        CovConfCases[[CvCase]]$muE <- Xmean + Xsd*StdDtRes[[CvCase]]$mu
        CovConfCases[[CvCase]]$ksiE <- Xmean + Xsd*StdDtRes[[CvCase]]$ksi
      }  else if (Type=="HomMxt") {
        displcmnt <- -Xmean
        sclfactor <- 1./Xsd
        CovConfCases[[CvCase]]$muE <- scale( scale(StdDtRes[[CvCase]]$mu,center=FALSE,scale=sclfactor),center=displcmnt,scale=FALSE )
        CovConfCases[[CvCase]]$ksiE <- 
          scale( scale(StdDtRes[[CvCase]]$ksi,center=FALSE,scale=sclfactor),center=displcmnt,scale=FALSE )
        attr(CovConfCases[[CvCase]]$muE,"scaled:center") <- attr(CovConfCases[[CvCase]]$ksiE,"scaled:center") <-
          attr(CovConfCases[[CvCase]]$muE,"scaled:scale") <- attr(CovConfCases[[CvCase]]$ksiE,"scaled:scale") <- NULL
      }  
      CovConfCases[[CvCase]]$SigmaE <- outer(Xsd,Xsd)*StdDtRes[[CvCase]]$Sigma
      CovConfCases[[CvCase]]$gamma1E <- StdDtRes[[CvCase]]$gamma1
      CovConfCases[[CvCase]]$OmegaE <- outer(Xsd,Xsd)*StdDtRes[[CvCase]]$Omega
      CovConfCases[[CvCase]]$alphaE <- StdDtRes[[CvCase]]$alpha
      logLiks[CvCase] <- CovConfCases[[CvCase]]$logLik <- StdDtRes[[CvCase]]$lnLik + lglikdif

      if (Type=="SingDst")  {
        names(CovConfCases[[CvCase]]$muE) <-  names(CovConfCases[[CvCase]]$ksiE) <- Xnames
      }  else { 
        colnames(CovConfCases[[CvCase]]$muE) <-  colnames(CovConfCases[[CvCase]]$ksiE) <- Xnames
        rownames(CovConfCases[[CvCase]]$muE) <-  rownames(CovConfCases[[CvCase]]$ksiE) <- levels(grouping)
        rownames(CovConfCases[[CvCase]]$muE) <-  rownames(CovConfCases[[CvCase]]$ksiE) <- grplvls
      }
      names(CovConfCases[[CvCase]]$gamma1E) <- names(CovConfCases[[CvCase]]$alphaE) <-
        dimnames(CovConfCases[[CvCase]]$SigmaE)[[1]] <- dimnames(CovConfCases[[CvCase]]$SigmaE)[[2]] <- 
        dimnames(CovConfCases[[CvCase]]$OmegaE)[[1]] <- dimnames(CovConfCases[[CvCase]]$OmegaE)[[2]] <- Xnames
      AICs[CvCase] <- CovConfCases[[CvCase]]$AIC <- -2*CovConfCases[[CvCase]]$logLik + 2*SKnpar(Conf,p,q,Ngrps=k)
      BICs[CvCase] <- CovConfCases[[CvCase]]$BIC <- -2*CovConfCases[[CvCase]]$logLik + log(n)*SKnpar(Conf,p,q,Ngrps=k)
    }
    if ( (StdDtRes[[CvCase]]$c2 > bordertol) || (maxsk-max(abs(StdDtRes[[CvCase]]$gamma1)) < bordertol) )
    {
      vcovl <- Getvcov(StdDtRes[[CvCase]])
      CovConfCases[[CvCase]]$status <- vcovl$status
      CovConfCases[[CvCase]]$mleCPvcov <- vcovl$mleCPvcov
      CovConfCases[[CvCase]]$muEse <- vcovl$muEse
      CovConfCases[[CvCase]]$SigmaEse <- vcovl$SigmaEse
      CovConfCases[[CvCase]]$gamma1Ese <- vcovl$gamma1Ese
    } else {
      CovConfCases[[CvCase]]$status <- "Onborder"
    }
  }
  if (SelCrit=="AIC")  {
    bestmod <- which.min(AICs)
  }  else if (SelCrit=="BIC") {
    bestmod <- which.min(BICs)
  }

  if (Type=="SingDst")  {
    return( new("IdtSngSNDE",ModelNames=names(AICs),ModelType=rep("SkewNormal",nCovCases),ModelConfig=1:nCovCases,
      NIVar=q,SelCrit=SelCrit,CovConfCases=CovConfCases,logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=TRUE) )
  }
  if (Type=="HomMxt")  {
    return( new("IdtMxSNDE",Hmcdt=TRUE,ModelNames=names(AICs),ModelType=rep("SkewNormal",nCovCases),ModelConfig=1:nCovCases,
      NIVar=q,SelCrit=SelCrit,CovConfCases=CovConfCases,logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=k) )
  }
}

IdtFDMxtSNmle <- function(Idt, grouping, CVtol=1.0e-5, OptCntrl=list(), CovCaseArg, Config, SelCrit)
{
  n <- Idt@NObs
  q <- Idt@NIVar
  p <- 2*q
  nk <- as.numeric(table(grouping))
  k <- length(nk) 
  if (k==1)  {
    stop("The data belongs to one single group. A partition into at least two different groups is required\n")
  }
  Xnams <- c(names(Idt@MidP),names(Idt@LogR))
  lev <- levels(grouping)
  anams <- list(Xnams,Xnams,lev)
  mnams <- list(lev,Xnams)
  if (CovCaseArg)  {
    nCovCases <- 4 
    CovCaseMap <- c(1,NA,2,3,4)
    slotnames <- paste("SNModCovC",1:4,sep="")
  } else {
    nCovCases <- 5 
    CovCaseMap <- 1:5
    slotnames <- paste("SNC",1:5,sep="")
  }
  logLiks <- AICs <- BICs <- rep(NA_real_,nCovCases)
  CovConfCases <- vector("list",nCovCases)
  names(logLiks) <- names(AICs) <- names(BICs) <- names(CovConfCases) <- slotnames
  for (model in Config)
  {
    nparbyg <- SKnpar(model,p,q)
    CvCase <- CovCaseMap[model]	
    CovConfCases[[CvCase]] <- list(
      muE=matrix(nrow=k,ncol=p,dimnames=mnams),muEse=matrix(nrow=k,ncol=p,dimnames=mnams),
      ksiE=matrix(nrow=k,ncol=p,dimnames=mnams),
      gamma1E=matrix(nrow=k,ncol=p,dimnames=mnams),gamma1Ese=matrix(nrow=k,ncol=p,dimnames=mnams),
      alphaE=matrix(nrow=k,ncol=p,dimnames=mnams),
      SigmaE=array(dim=c(p,p,k),dimnames=anams),SigmaEse=array(dim=c(p,p,k),dimnames=anams),
      OmegaE=array(dim=c(p,p,k),dimnames=anams),
      mleCPvcov=array(dim=c(nparbyg,nparbyg,k),dimnames=list(NULL,NULL,levels(grouping))),
      status=NULL,logLik=0.,AIC=NULL,BIC=NULL,optres=vector("list",k)
    )
    names(CovConfCases[[CvCase]]$optres) <- lev
  }
  for (g in 1:k) {
    Idtg <- Idt[grouping==lev[g],]
    IdtgDF <- cbind(Idtg@MidP,Idtg@LogR)
    Xbar <- colMeans(IdtgDF)
    Xstdev <- sapply(IdtgDF,sd)
    CnstV <- which(Xstdev/abs(Xbar)<CVtol)
    if (length(CnstV)==1)  { 
      stop("Variable ",names(CnstV)," appears to be constant in group",lev[g],"\n")
    }  else if (length(CnstV)>1)  {
      stop("Variables ",paste(names(CnstV),collapse=" ")," appear to be constant in group ",lev[g],"\n")
    }
    pres <- IdtSNmle(Idtg,Type="SingDst",CVtol=CVtol,OptCntrl=OptCntrl,
      CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    for (model in Config)
    { 
      CvCase <- CovCaseMap[model]	
      CovConfCases[[CvCase]]$muE[g,] <- pres@CovConfCases[[CvCase]]$muE
      CovConfCases[[CvCase]]$ksiE[g,] <- pres@CovConfCases[[CvCase]]$ksiE
      CovConfCases[[CvCase]]$gamma1E[g,] <- pres@CovConfCases[[CvCase]]$gamma1E
      CovConfCases[[CvCase]]$alphaE[g,] <- pres@CovConfCases[[CvCase]]$alphaE
      CovConfCases[[CvCase]]$SigmaE[,,g] <- pres@CovConfCases[[CvCase]]$SigmaE
      CovConfCases[[CvCase]]$OmegaE[,,g] <- pres@CovConfCases[[CvCase]]$OmegaE
      CovConfCases[[CvCase]]$logLik <- CovConfCases[[CvCase]]$logLik + pres@CovConfCases[[CvCase]]$logLik
      CovConfCases[[CvCase]]$optres[[g]] <- pres@CovConfCases[[CvCase]]$optres
      if (!is.null(pres@CovConfCases[[CvCase]]$muEse))
        { CovConfCases[[CvCase]]$muEse[g,] <- pres@CovConfCases[[CvCase]]$muEse }
      if (!is.null(pres@CovConfCases[[CvCase]]$gamma1Ese))
        { CovConfCases[[CvCase]]$gamma1Ese[g,] <- pres@CovConfCases[[CvCase]]$gamma1Ese }
      if (!is.null(pres@CovConfCases[[CvCase]]$SigmaEse))
        { CovConfCases[[CvCase]]$SigmaEse[,,g] <- pres@CovConfCases[[CvCase]]$SigmaEse }
      if (!is.null(pres@CovConfCases[[CvCase]]$mleCPvcov))  {
        CovConfCases[[CvCase]]$mleCPvcov[,,g] <- pres@CovConfCases[[CvCase]]$mleCPvcov
        if (g==1) {
          dimnames(CovConfCases[[CvCase]]$mleCPvcov)[[1]]  <- dimnames(CovConfCases[[CvCase]]$mleCPvcov)[[2]] <- 
            rownames(pres@CovConfCases[[CvCase]]$mleCPvcov)
        }
      }
    }            
  }
  for (model in Config)
  {
    CvCase <- CovCaseMap[model]	
    logLiks[CvCase] <- CovConfCases[[CvCase]]$logLik 
    nmodelfreepar <- SKnpar(model,p,q,Ngrps=k,Mxt="GenMod")
    AICs[CvCase] <- CovConfCases[[CvCase]]$AIC <- -2*CovConfCases[[CvCase]]$logLik + 2 * nmodelfreepar
    BICs[CvCase] <- CovConfCases[[CvCase]]$BIC <- -2*CovConfCases[[CvCase]]$logLik + log(n) * nmodelfreepar
  }
  if (SelCrit=="AIC")  {
    bestmod = which.min(AICs)
  }  else if (SelCrit=="BIC")  {
    bestmod = which.min(BICs)
  }
 
  new("IdtMxSNDE",ModelNames=names(AICs),ModelType=rep("SkewNormal",nCovCases),ModelConfig=1:nCovCases,
    grouping=grouping,Hmcdt=FALSE,CovConfCases=CovConfCases,SelCrit=SelCrit,NIVar=q,
    logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=k)
}

