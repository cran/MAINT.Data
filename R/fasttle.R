setMethod("fasttle",
  signature(Idt = "IData"),
  function(Idt,
    CovCase=1:4,
    SelCrit=c("BIC","AIC"),
    alpha=control@alpha,
    nsamp = control@nsamp,
    seed=control@seed,
    trace=control@trace,
    use.correction=control@use.correction,
    ncsteps=control@ncsteps,
    getalpha=control@getalpha,
    rawMD2Dist=control@rawMD2Dist,				
    MD2Dist=control@MD2Dist,
    eta=control@eta,
    multiCmpCor=control@multiCmpCor,				
    getkdblstar=control@getkdblstar,
    outlin=control@outlin,
    trialmethod=control@trialmethod,
    m=control@m,
    reweighted = control@reweighted,
    otpType=control@otpType,
    control=RobEstControl(), ...)
  {
    SelCrit <- match.arg(SelCrit)
    if (!requireNamespace("robustbase",quietly=TRUE)) {
      stop("fasttle needs the robustbase package to work. Please install it.\n")
    }
    if  (getkdblstar=="Twopplusone") { 
      kdblstar <- 2*Idt@NIVar+1
    }  else {
      if (!is.finite(getkdblstar)) {
        stop("Wrong value for argument getkdblstar\n")
      }
      kdblstar <- getkdblstar 
    }
    q <- Idt@NIVar
    if (getalpha=="TwoStep")
    {
      X <- cbind(Idt@MidP,Idt@LogR)
      if (MD2Dist=="ChiSq") {
        fstsol <- fasttle1(Idt,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
          rawMD2Dist,eta,multiCmpCor,kdblstar,outlin,trialmethod,m,reweighted,otpType="OnlyEst")
        nOtls <- MDOtlDet(X,coef(fstsol)$mu,coef(fstsol)$Sigma,eta=eta,RefDist="ChiSq",multiCmpCor=multiCmpCor,otp="onlycnt")
      }  else if (MD2Dist=="CerioliBetaF") {
        fstsol <- fasttle1(Idt,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
          rawMD2Dist,eta,multiCmpCor,kdblstar,outlin,trialmethod,m,reweighted,otpType="SetMD2andEst")
        nOtls <- MDOtlDet(X,coef(fstsol$sol)$mu,coef(fstsol$sol)$Sigma,eta=eta,RefDist="CerioliBetaF",
          Rewind=fstsol$RewghtdSet,multiCmpCor=multiCmpCor,otp="onlycnt")
      }
      if (nOtls==0) {
        Config <- getConfig(...)
        if (is.null(Config))  
        {
          Config <- ifelse(CovCase==1,1,CovCase+1)
          CovCaseArg <- TRUE	
        } else {  
          CovCaseArg <- FALSE
        }	
        finalsol <- IdtNmle(Idt,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
        warning(paste("fasttle returned the classical maximum likelihood estimates because the data does not appear to include any outlier.\n",
          "If you want to force a trimmed likelihood estimator run fasttle with the argument getalpha=FALSE.\n"))
        if (otpType=="OnlyEst") {
          return(finalsol)
        }  else {
          rawSet <- RewghtdSet <- 1:Idt@NObs
          names(rawSet) <- names(RewghtdSet) <- Idt@ObsNames
          if (outlin=="MidPandLogR") {
            RobMD2 <- GetMD2(X,coef(finalsol)$mu,coef(finalsol)$Sigma)
          } else if (outlin=="MidP") {
            RobMD2 <- GetMD2(X[,1:q],coef(finalsol)$mu[1:q],coef(finalsol)$Sigma[1:q,1:q])
          } else if (outlin=="LogR") {
            RobMD2 <- GetMD2(X[,(q+1):(2*q)],coef(finalsol)$mu[(q+1):(2*q)],coef(finalsol)$Sigma[(q+1):(2*q),(q+1):(2*q)])
          } 
          if (otpType=="SetMD2andEst") { 
            return(list(sol=finalsol,rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=RobMD2,
              cnp2=c(1.,1.),raw.cov=coef(fstsol)$Sigma,raw.cnp2=c(1.,1.)))
          }  else if (otpType=="SetMD2EstandPrfSt") {
            return(list(sol=finalsol,rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=RobMD2,
              cnp2=c(1.,1.),raw.cov=coef(fstsol)$Sigma,raw.cnp2=c(1.,1.),
              PerfSt=list(RepSteps=NULL,RepLogLik=NULL,StpLogLik=NULL)))
          }  
        }
      }  
      newalpha <- 1. - nOtls/Idt@NObs
      return( fasttle1(Idt,CovCase,SelCrit,newalpha,nsamp,ncsteps,trace,use.correction,
        rawMD2Dist,eta,multiCmpCor,kdblstar,outlin,trialmethod,m,reweighted,otpType) )
    }  else {
      return( fasttle1(Idt,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
        rawMD2Dist,eta,multiCmpCor,kdblstar,outlin,trialmethod,m,reweighted,otpType) )
    } 
  }
)

RobEstControl <- function (alpha=0.75,
                           nsamp=500,
                           seed=NULL,
                           trace=FALSE,
                           use.correction=TRUE,
                           ncsteps=200,
                           getalpha = "TwoStep",
                           rawMD2Dist="ChiSq",				
                           MD2Dist="ChiSq",
                           eta=0.025,   
                           multiCmpCor="never",				
                           getkdblstar="Twopplusone",
                           outlin="MidPandLogR",
                           trialmethod="simple",
                           m=1,
                           reweighted=TRUE,
                           otpType="OnlyEst")
{
    new("RobEstControl", alpha = alpha,
                         nsamp = nsamp,
                         seed = seed,
                         trace = trace,
                         use.correction = use.correction,
                         ncsteps=200,
                         getalpha = getalpha,
                         rawMD2Dist=rawMD2Dist,				
                         MD2Dist=MD2Dist,				
                         eta=eta,   
                         multiCmpCor=multiCmpCor,				
                         getkdblstar = getkdblstar,
                         outlin =outlin,
                         trialmethod = trialmethod,
                         m=m,
                         reweighted=reweighted,
                         otpType = otpType )
}

fasttle1 <- function(data,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
  rawMD2Dist,eta,multiCmpCor,kdblstar,outlin,trialmethod,m,reweighted,otpType,Vnames=NULL,...)
{
  datatype <- class(data)
  if (datatype!="IData" && datatype!="matrix" && datatype!="data.frame")  {
    stop("Wrong class for data argument\n")
  }

  Config <- getConfig(...)
  if (is.null(Config))  
  {
    Config <- ifelse(CovCase==1,1,CovCase+1)
    CovCaseArg <- TRUE	
    nCovCases <- 4 
    CovCaseMap <- c(1,NA,2,3,4)
    modnames <- paste("NModCovC",1:4,sep="")
  } else {  
    CovCaseArg <- FALSE
    nCovCases <- 5 
    CovCaseMap <- 1:5
    modnames <- paste("NC",1:5,sep="")
  }	
  bestSetbyCvC <- vector("list",nCovCases)
  bestCrt <- rep(Inf,nCovCases)
  if (trialmethod=="Poolm") {
    Poolm <- 1
  }  else {
    Poolm <- 0
  }
  if ((Poolm==1) && m==1) {
      stop("m argument needs to be set to a value higher than 1 when the 'Poolm' trialmethod is chosen.\n") 
  }
  if (otpType=="SetMD2EstandPrfSt") {
    ClctSt <- 1
    if (outlin=="MidPandLogR")  {
      RepSteps <- vector("list",nCovCases)
      RepLogLik <- vector("list",nCovCases)
      StpLogLik <- vector("list",nCovCases)
      names(RepSteps) <- names(RepLogLik) <- names(StpLogLik) <- modnames
    }
    Repnames <- paste("Rep",1:nsamp,sep="")
    Stepnames <- paste("Stp",1:ncsteps,sep="")
  }  else {
    ClctSt <- 0
  }

  if (datatype=="IData")
  {
    n <- data@NObs  
    if (outlin=="MidPandLogR")  {
      X <- as.matrix(cbind(data@MidP,data@LogR))
      p <- 2*data@NIVar
    }  else {
      p <- data@NIVar
      Config1 <- unique(ifelse(is.element(Config,c(1,4)),1,5)) 
      if (outlin=="MidP") {
        X <- as.matrix(data@MidP)
        Vind <- 1:p
      }
      if (outlin=="LogR") {
        X <- as.matrix(data@LogR)
        Vind <- (p+1):(2*p)
      }
    }
    rownames(X) <- data@ObsNames
  }  else  {
    X <- as.matrix(data)
    n <- nrow(X)
    p <- ncol(X)
  }  

  k <- robustbase::h.alpha.n(alpha,n,p)
    
  c0 <- -0.5*(p*log(2*pi))
  if (SelCrit=="BIC") {
    penC <- log(n)
  }  else if (SelCrit=="AIC") {
    penC <- 2
  }
  if (outlin=="MidPandLogR")  {
    for (Cnf in Config) {
      CvCase <- CovCaseMap[Cnf]	
      if (Cnf==2)  { 
        if (datatype!="IData")  {
          stop("Wrong class for data argument\n")
        }
        Cftmpsol <- Rfasttle(data,kdblstar,k,nsamp,2,SelCrit,maxrefstps=ncsteps,...)
      }  else {
        Cftmpsol <- .Call( "Cfasttle", X, n, p, Poolm, m, kdblstar, k, nsamp, Cnf,
          c0, ncsteps, ClctSt, PACKAGE = "MAINT.Data" )
        Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
        if (datatype=="IData") {
          names(Cftmpsol$Set) <- data@ObsNames[Cftmpsol$Set]
        } else {
          names(Cftmpsol$Set) <- rownames(data)[Cftmpsol$Set]
        }
      }
      CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(Cnf,p,p/2)
      if (CmpCrt < bestCrt[CvCase]) {
        bestCrt[CvCase] <- CmpCrt
        bestSetbyCvC[[CvCase]] <- Cftmpsol$Set
      }
      if (otpType=="SetMD2EstandPrfSt") {
        RepSteps[[CvCase]] <- Cftmpsol$RepSteps
        maxnSteps <- max(RepSteps[[CvCase]])
        RepLogLik[[CvCase]] <- Cftmpsol$RepLogLik
        names(RepSteps[[CvCase]]) <- names(RepLogLik[[CvCase]]) <- Repnames
        StpLogLik[[CvCase]] <- Cftmpsol$StpLogLik[,1:maxnSteps,drop=FALSE]
        dimnames(StpLogLik[[CvCase]]) <- list(Repnames,Stepnames[1:maxnSteps])
      }
    } 
  } else {
    for (Cnf in Config1)
    {
      CvCase <- CovCaseMap[Cnf]	
      Cftmpsol <- .Call( "Cfasttle", X, n, p, Poolm, m, kdblstar, k, nsamp, CvCase,
        c0, ncsteps, ClctSt, PACKAGE = "MAINT.Data" )
      Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
      CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(Cnf,2*p,p)
      if (otpType=="SetMD2EstandPrfSt") {
        RepSteps[[CvCase]] <- Cftmpsol$RepSteps
        maxnSteps <- max(RepSteps[[CvCase]])
        RepLogLik[[CvCase]] <- Cftmpsol$RepLogLik
        names(RepSteps[[CvCase]]) <- names(RepLogLik[[CvCase]]) <- Repnames
        StpLogLik[[CvCase]] <- Cftmpsol$StpLogLik[,1:maxnSteps]
        dimnames(StpLogLik[[CvCase]]) <- list(Repnames,Stepnames[1:maxnSteps])
      }
      if (CmpCrt < bestCrt[CvCase]) {
        bestCrt[CvCase] <- CmpCrt
        bestSetbyCvC[[CvCase]] <- Cftmpsol$Set
      }
    }
  }

  bestSet <- bestSetbyCvC[[which.min(bestCrt)]]
  if (datatype!="IData") { data <- IData(as.data.frame(X),"AllMidP_AllLogR",VarNames=Vnames) }
  bestsol <- IdtNmle(data[bestSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
  
  dhn1 <- robustbase::.MCDcons(p,k/n)
  if (use.correction) {
    dhn <- dhn1*MCDcnp2(p,n,alpha)
  } else {
    dhn <- rep(dhn1,4)
  }
  if (!reweighted)  {
    if (outlin=="MidPandLogR") {
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$mleSigE <- dhn[CvCase] * bestsol@CovConfCases[[CvCase]]$mleSigE
      }
    } else {
      for (Cnf in Config) {
         CvCase <- CovCaseMap[Cnf]	
         bestsol@CovConfCases[[CvCase]]$mleSigE[Vind,Vind] <- dhn[CvCase] * bestsol@CovConfCases[[CvCase]]$mleSigE[Vind,Vind]
      }
    }  
    RewghtdSet <- NULL
  } else {
    corrawcov <- bestsol@CovConfCases[[bestsol@BestModel]]$mleSigE
    if (outlin=="MidPandLogR")
    {
      Xdev <- scale(X,center=bestsol@mleNmuE,scale=FALSE)
      Sigma <- corrawcov <- dhn[bestsol@BestModel] * corrawcov       
    } else { 
      Xdev <- scale(X,center=bestsol@mleNmuE[Vind],scale=FALSE)
      Sigma <- corrawcov[Vind,Vind] <- dhn[bestsol@BestModel] * corrawcov[Vind,Vind]
    }
    SigmaI <- pdwt.solve(Sigma)
    MD2 <- apply(Xdev,1,function(x) x%*%SigmaI%*%x)

    oneminuseta <- 1-eta
    if (multiCmpCor=="always" || multiCmpCor=="iterstep") {
      oneminusalpha <- oneminuseta^(1/n)
    } else {
      oneminusalpha <- oneminuseta
    }
    h <- length(bestSet)
    findotl <- TRUE
    iter <- 1
    while (findotl && iter<3)  { 
      if (rawMD2Dist=="ChiSq")  { 
        MD2trshld <- qchisq(oneminusalpha,p)
      } else if (rawMD2Dist=="HardRockeAdjF")  {
        MD2trshld <-qHardRoqF(oneminusalpha,n,p,h) / dhn[bestsol@BestModel]
      } else if (rawMD2Dist=="HardRockeAsF")  {
        MD2trshld <- qHardRoqF(oneminusalpha,n,p,h,adj=FALSE) / dhn[bestsol@BestModel]
      }
      RewghtdSet <- which(MD2<=MD2trshld)
      if (multiCmpCor=="iterstep")
      {
        if (length(RewghtdSet)==n) {
          findotl <- FALSE
        } else {
          oneminusalpha <- oneminuseta
          iter <- iter+1
        }
      } else {
        findotl <- FALSE
      }
    }

    bestsol <- IdtNmle(data[RewghtdSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    m <- length(RewghtdSet)
    rdmhn1 <- robustbase::.MCDcons(p,m/n)
    if (outlin=="MidPandLogR")  {
      if (use.correction) {
        rdmhn <- rdmhn1 * MCDRewcnp2(p,n,m/n,alpha)
      } else { 
        rdmhn <- rep(rdmhn1,4)
      }
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$mleSigE <- rdmhn[CvCase] * bestsol@CovConfCases[[CvCase]]$mleSigE
      }
    } else {
      if (use.correction) {
        rdmhn <- rdmhn1 * MCDRewcnp2(data@NIVar,n,m/n,alpha)
      } else { 
        rdmhn <- rep(rdmhn1,4)
      }
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$mleSigE[Vind,Vind] <- rdmhn[CvCase] * bestsol@CovConfCases[[CvCase]]$mleSigE[Vind,Vind]
      }
    }
  }

  finalsol <- new("IdtSngNDRE",ModelNames=bestsol@ModelNames,ModelType=bestsol@ModelType,ModelConfig=bestsol@ModelConfig,
    NIVar=bestsol@NIVar,SelCrit=bestsol@SelCrit,logLiks=bestsol@logLiks,BICs=bestsol@BICs,AICs=bestsol@AICs,
    BestModel=bestsol@BestModel,RobNmuE=bestsol@mleNmuE,CovConfCases=bestsol@CovConfCases,SngD=TRUE)
  for (case in 1:length(finalsol@CovConfCases)) {
    if (!is.null(finalsol@CovConfCases[[case]])) {
      names(finalsol@CovConfCases[[case]])[1] <- "RobSigE"
      finalsol@CovConfCases[[case]][2] <- finalsol@CovConfCases[[case]][3] <- NULL
    }
  }
  
  if (otpType=="OnlyEst") {
    return(finalsol)
  }  else {
    rawSet <- sort(bestSet)
    cnp2 <- c(rdmhn1,rdmhn[bestsol@BestModel]/rdmhn1)
    raw.cnp2 <- c(dhn1,dhn[bestsol@BestModel]/dhn1)
    names(cnp2) <- names(raw.cnp2) <- NULL 
    if (otpType=="SetMD2andEst") { 
      return(list(sol=finalsol,rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=MD2,
        cnp2=cnp2,raw.cov=corrawcov,raw.cnp2=raw.cnp2))
    }  else if (otpType=="SetMD2EstandPrfSt") {
      return(list(sol=finalsol,rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=MD2,
        cnp2=cnp2,raw.cov=corrawcov,raw.cnp2=raw.cnp2,
        PerfSt=list(RepSteps=RepSteps,RepLogLik=RepLogLik,StpLogLik=StpLogLik)))
    }
  }
}

getIdtOutl <- function(Idt,IdtE=NULL,muE=NULL,SigE=NULL,
  eta=0.025,Rewind=IdtE$RewghtdSet,m=length(Rewind),
  RefDist=c("ChiSq","HardRockeAdjF","HardRockeAsF","CerioliBetaF"),
  multiCmpCor=c("never","always","iterstep"),outlin=c("MidPandLogR","MidP","LogR"))
{
  RefDist <- match.arg(RefDist) 
  multiCmpCor <- match.arg(multiCmpCor)
  outlin <- match.arg(outlin)

  if ( is.null(muE) || is.null(SigE) ) {
    if (is.null(IdtE)) {
      stop("Missing mean and/or covariance estimates in call to getIdtOutl.\n")
    }  
    if (!isClass(IdtE$sol,"IdtSngNDE")) {
      stop("sol component of IdtE argument is not an object of class IdtSngNDRE or IdtSngNDE as required.\n")
    }
    if (is.null(muE)) muE <- coef(IdtE$sol)$mu
    if (is.null(SigE)) SigE <- coef(IdtE$sol)$Sigma
  }

  X <- data.frame(cbind(Idt@MidP,Idt@LogR))
  if (outlin=="MidPandLogR") vind <- 1:(2*Idt@NIVar)
  else if (outlin=="MidP") vind <- 1:Idt@NIVar
  else if (outlin=="LogR") vind <- (Idt@NIVar+1):(2*Idt@NIVar)

  MDOtlDet(X[vind],muE[vind],SigE[vind,vind],eta,m,ret="Outliers",RefDist=RefDist,Rewind=Rewind,multiCmpCor=multiCmpCor)
}

GetMD2 <- function(Data,muE,SigE)
{
  SigInv <- solve(SigE)
  XDev <- scale(Data,center=muE,scale=FALSE)
  apply(XDev,1,function(x) x %*% SigInv %*% x)
}

MDOtlDet <- function(Data,muE,SigE,eta,h=NULL,ret=c("Outliers","Regular"),
  RefDist=c("ChiSq","HardRockeAdjF","HardRockeAsF","CerioliBetaF"),
  Rewind=NULL,multiCmpCor=c("never","always","iterstep"),otp=c("indices","onlycnt","indandMD2"))
{
  ret <- match.arg(ret)
  otp <- match.arg(otp)
  RefDist <- match.arg(RefDist)
  multiCmpCor <- match.arg(multiCmpCor)
  
  if ( (RefDist=="HardRockeAdjF" || RefDist=="HardRockeAdjF") && is.null(h) ) {  
    stop("In order to use the Hardin and Rocke F approximations you have to specify the number\n",
      "of the observations kept in the trimmed raw MCD estimator with the h argument.") 
  }
  if (RefDist=="CerioliBetaF") {
    if (is.null(Rewind)) {  
      stop("In order to use the Cerioli Beta and F approximation you have to specify the indices\n",
        "of the observations kept in the reweighted step using the Rewind argument.")
    }
    h <- length(Rewind)
  } 
  if (!is.matrix(Data)) Data <- as.matrix(Data)
  n <- nrow(Data)
  p <- ncol(Data)
  if (length(muE)!=p) stop("Wrong muE dimension\n")
  if (!is.matrix(SigE)) SigE <- as.matrix(SigE)
  if (nrow(SigE)!=p || ncol(SigE)!=p) stop("Wrong SigE dimension\n")
	
  MDist <- GetMD2(Data,muE,SigE)
  oneminuseta <- 1-eta
  if (multiCmpCor=="always" || multiCmpCor=="iterstep") {
    oneminusalpha <- oneminuseta^(1/n)
  } else {
    oneminusalpha <- oneminuseta
  }

  findotl <- TRUE
  iter <- 1
  while (findotl && iter<3)  { 
    if (RefDist=="ChiSq")  { 
      delta <- qchisq(oneminusalpha,p)
    } else if (RefDist=="HardRockeAdjF")  {
      delta <- qHardRoqF(oneminusalpha,n,p,h)
    } else if (RefDist=="HardRockeAsF")  {
      delta <- qHardRoqF(oneminusalpha,n,p,h,adj=FALSE)
    } else if (RefDist=="CerioliBetaF")  {
      delta1 <- ((h-1)^2/h) * qbeta(oneminusalpha,p/2,(h-p-1)/2)
      delta2 <- (((h+1)*(h-1)*p)/(h*(h-p))) * qf(oneminusalpha,p,h-p)
    } 
    if (ret=="Outliers")  {
      if (RefDist!="CerioliBetaF") {
        otl <- which(MDist>delta)
      } else {
        boolRewind <- sapply(1:n,function(x) is.element(x,Rewind))
        otl <- which(ifelse(boolRewind,MDist>delta1,MDist>delta2))
      }  
      if (length(otl)==0) {
        otl <- NULL
      }
      if (multiCmpCor=="iterstep") {
        if (length(otl)==0) {
          findotl <- FALSE
        } else {
          oneminusalpha <- oneminuseta
        }
      }
    }  else if (ret=="Regular") {
      if (RefDist!="CerioliBetaF") {
        reg <- which(MDist<=delta)
      } else {
        boolRewind <- sapply(1:n,function(x) is.element(x,Rewind))
        reg <- which(ifelse(boolRewind,MDist<=delta1,MDist<=delta2))
      }  
      if (length(reg)==0) {
        reg <- NULL
      }
      if (multiCmpCor=="iterstep") {
        if (length(reg)==n) {
          findotl <- FALSE
        } else {
          oneminusalpha <- oneminuseta
        }
      }
    }
    if (multiCmpCor!="iterstep") findotl <- FALSE
    iter <- iter+1
  }

  if (ret=="Outliers")  {
    if (otp=="indices")  {
        return(otl)
      } else if (otp=="onlycnt")  {
        return(length(otl))
      } else if (otp=="indandMD2")  {
        return( list(outliers=otl,MD2=MDist) )
      }      
  } else if (ret=="Regular") {
    if (otp=="indices")  {
      return(reg)
    } else if (otp=="onlycnt")  {
      return(length(reg))
    } else if (otp=="indandMD2")  {
      return( list(regobs=reg,MD2=MDist) )
    }
  }
}

MCDcnp2 <- function(p,n,alpha,userobustbase=FALSE)
{
  dhn2i <- numeric(2)
  lnp <- log(p)
  A <- matrix(c(1.,1.,-log(3*p^2),-log(5*p^2)),2,2)

  MCDcnp21 <- function(Cfind,p,n,w)
  {
    if (userobustbase && Cfind==1) return(robustbase::.MCDcnp2(p,n,alpha))
    for (alphaind in 1:2)
    {	
      x <- solve(A,log(-SigEf_etas[,alphaind,Cfind])-SigEf_kappas[,alphaind,Cfind]*lnp)
      gambetest <- c(-exp(x[1]),x[2])
      dhn2i[alphaind] <- 1. + gambetest[1]/n^gambetest[2]
    }
    1./ ( (1.-w)*dhn2i[1]+w*dhn2i[2] )
  }  

  sapply( 1:4,MCDcnp21,p=p,n=n, w=(alpha-SigEf_alphas[1])/(SigEf_alphas[2]-SigEf_alphas[1]) )
}


MCDRewcnp2 <- function(p,n,movern,alpha,userobustbase=FALSE)
{ 
  Refp2alphahat <- function(Cfind)
  {
    if (userobustbase && Cfind==1) return(1./robustbase::.MCDcnp2.rew(p,n,alpha))
    uplims[Cfind] + SigRewEf_gammas[1,Cfind]/n^SigRewEf_betas[1,Cfind]
  }

  MCDRewcnp21 <- function(Cfind)
  {
    if (userobustbase && Cfind==1) return(1./robustbase::.MCDcnp2.rew(p,n,alpha))
    x <- solve(A,log(-SigRewEf_etas[,Cfind])-SigRewEf_kappas[,Cfind]*lnp)
    gambetest <- c(-exp(x[1]),x[2])
    uplims[Cfind] + gambetest[1]/n^gambetest[2]
  }  
  
  uplims <- 1. - SigRewEf_Betams - SigRewEf_Betaalphas
  if (p==2) Betas0 <- sapply(1:4,Refp2alphahat)
  else {
    lnp <- log(p)
    A <- matrix(c(1.,1.,-log(3*p^2),-log(5*p^2)),2,2)
    Betas0 <- sapply(1:4,MCDRewcnp21)
  }  
  if (userobustbase) SigRewEf_Betams[1] <- SigRewEf_Betaalphas[1] <- 0.  

  1. / (Betas0 + SigRewEf_Betams*movern + SigRewEf_Betaalphas*alpha)
}
