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
      stop("fasttle needs the robustbase package to work. Please install it\n")
    }
    if  (getkdblstar=="Twopplusone") { 
      kdblstar <- 2*Idt@NIVar+1
    }  else {
      if (!is.finite(getkdblstar)) {
        stop("Wrong value for argument getkdblstar\n")
      }
      kdblstar <- getkdblstar 
    }
    if (getalpha=="TwoStep")
    {
      fstsol <- fasttle1(Idt,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
        kdblstar,outlin,trialmethod,m,reweighted,otpType="OnlyEst")
      nOtls <- MDOtlDet(data.frame(cbind(Idt@MidP,Idt@LogR)),coef(fstsol)$mu,coef(fstsol)$Sigma,0.025,onlycnt=TRUE) 
      newalpha <- 1. - nOtls/Idt@NObs
      return( fasttle1(Idt,CovCase,SelCrit,newalpha,nsamp,ncsteps,trace,use.correction,
        kdblstar,outlin,trialmethod,m,reweighted,otpType) )
    }  else {
      return( fasttle1(Idt,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
        kdblstar,outlin,trialmethod,m,reweighted,otpType) )
    } 
  }
)

MDOtlDet <- function(Data,muE,SigE,alpha,ret=c("Otliers","Regular"),onlycnt=FALSE)
{
  ret <- match.arg(ret)
  if (!is.matrix(Data)) Data <- as.matrix(Data)
  n <- nrow(Data) ;  p <- ncol(Data)
  if (length(muE)!=p) stop("Wrong muE dimension\n")
  if (!is.matrix(SigE)) SigE <- as.matrix(SigE)
  if (nrow(SigE)!=p || ncol(SigE)!=p) stop("Wrong SigE dimension\n")
	
  SigInv <- solve(SigE)
  XDev <- scale(Data,center=muE,scale=FALSE)
  MDist <- apply(XDev,1,function(x) x %*% SigInv %*% x)  # Things to do:  see if fastR package improves speed in here !!

  oneminusalpha <- 1-alpha
  delta <- qchisq(oneminusalpha,p)
  if (ret=="Otliers" && !onlycnt)  {
    return( which(MDist>delta) )
  } else if (ret=="Otliers" && onlycnt)  {
    return( length(which(MDist>delta)) )
  }  else if (ret=="Regular" && !onlycnt)  {  
    return( which(MDist<=delta) )
  }  else if (ret=="Regular" && onlycnt)  {  
    return( length(which(MDist<=delta)) )
  }
}

RobEstControl <- function (alpha=0.75,
                           nsamp=500,
                           seed=NULL,
                           trace=FALSE,
                           use.correction=TRUE,
                           ncsteps=200,
                           getalpha = "TwoStep",
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
                         getkdblstar = getkdblstar,
                         outlin =outlin,
                         trialmethod = trialmethod,
                         m=m,
                         reweighted=reweighted,
                         otpType = otpType )
}

fasttle1 <-   function(data,CovCase,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,
  kdblstar,outlin,trialmethod,m,reweighted,otpType,Vnames=NULL,...)
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
        StpLogLik[[CvCase]] <- Cftmpsol$StpLogLik[,1:maxnSteps]
        dimnames(StpLogLik[[CvCase]]) <- list(Repnames,Stepnames[1:maxnSteps])
      }
    }
  }  else {
    Cftmpsol <- .Call( "Cfasttle", X, n, p, Poolm, m, kdblstar, k, nsamp, 1,
      c0, ncsteps, ClctSt, PACKAGE = "MAINT.Data" )
    Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
    CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(1,2*p,p)
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

  bestSet <- bestSetbyCvC[[which.min(bestCrt)]]
  if (datatype!="IData") { data <- IData(as.data.frame(X),"AllMidP_AllLogR",VarNames=Vnames) }
  bestsol <- IdtNmle(data[bestSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)

  if (use.correction) {
    dhn1 <- robustbase::.MCDcons(p,k/n)
  }  else  {
    dhn1 <- 1.
  }
  if (!reweighted && use.correction) {
    dhn2 <- MCDcnp2(p,n,alpha,Config)
    if (outlin=="MidPandLogR")
    {
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$mleSigE <- dhn1*dhn2[Cnf] * bestsol@CovConfCases[[CvCase]]$mleSigE
      }
    } else {
      bestsol@CovConfCases[[1]]$mleSigE[Vind,Vind] <- dhn1*dhn2[1] * bestsol@CovConfCases[[1]]$mleSigE[Vind,Vind]
    } 
  }
  if (reweighted)
  {
    if (use.correction)
    {
      dhn <- dhn1 * MCDcnp2(p,n,alpha,bestsol@BestModel)
    } else { 
      dhn <- dhn1
    }
    if (outlin=="MidPandLogR")
    {
      Xdev <- scale(X,center=bestsol@mleNmuE,scale=FALSE)
      Sigma <- dhn * bestsol@CovConfCases[[bestsol@BestModel]]$mleSigE      
    } else { 
      Xdev <- scale(X,center=bestsol@mleNmuE[Vind],scale=FALSE)
      Sigma <- dhn * bestsol@CovConfCases[[1]]$mleSigE[Vind,Vind]
    }
    SigmaI <- pdwt.solve(Sigma)
    MD2 <- apply(Xdev,1,function(x) x%*%SigmaI%*%x)
    MD2trshld <- qchisq(0.975,p)
    RewghtdSet <- which(MD2<=MD2trshld)
    bestsol <- IdtNmle(data[RewghtdSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
  } else {
    RewghtdSet <- NULL
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
    if (otpType=="SetMD2andEst") { 
      return(list(sol=finalsol,rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=MD2))
    }  else if (otpType=="SetMD2EstandPrfSt") {
      return(list(sol=finalsol,rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=MD2,
        PerfSt=list(RepSteps=RepSteps,RepLogLik=RepLogLik,StpLogLik=StpLogLik)))
    }
  }
}

MCDcnp2 <- function(p,n,alpha,Cf)
{
  MCDcnp21 <- function(Cfind,pind,n,w)
  {
    gammas <- SigEf_gammas[pind,,Cfind]
    betas <- SigEf_betas[pind,,Cfind]
    dhn2i <- 1. + gammas/n^betas
    1./ ( (1.-w)*dhn2i[1]+w*dhn2i[2] )
  }  

  pind <- which(p==SigEf_ps)
  if (is.null(pind))
    stop(paste("Small sample correction was not implemented yet for p =",p,".\n"))
  sapply(ifelse(Cf==1,Cf,Cf-1),MCDcnp21,pind=pind,n=n,
    w=(alpha-SigEf_alphas[1])/(SigEf_alphas[2]-SigEf_alphas[1]))
}
