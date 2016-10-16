setMethod("fulltle",
  signature(Idt = "IData"),
  function(Idt, alpha=0.75, reweighted=TRUE, CorrF=c("smallsmp","consistent","none"),
    outlin=c("MidPandLogR","MidP","LogR"),
    CovCase=1:4, SelCrit=c("BIC","AIC"), force=FALSE, otpType=c("OnlyEst","SetMD2andEst"), ... )
  {
    if (!requireNamespace("robustbase",quietly=TRUE)) 
      stop("fulltle needs the robustbase package to work. Please install it\n")

    SelCrit <- match.arg(SelCrit)
    otpType <- match.arg(otpType)
    CorrF <- match.arg(CorrF)
    outlin <- match.arg(outlin)

    Config <- getConfig(...)
    if (is.null(Config))  
    {
      Config <- ifelse(CovCase==1,1,CovCase+1)
      CovCaseArg <- TRUE	
    } else {  
      CovCaseArg <- FALSE
    }	

    n <- Idt@NObs    
    if (outlin=="MidPandLogR")  {
      X <- as.matrix(cbind(Idt@MidP,Idt@LogR))
      p <- 2*Idt@NIVar
    }
    else {
      p <- Idt@NIVar
      if (outlin=="MidP") {
        X <- as.matrix(Idt@MidP)
        Vind <- 1:p
      }
      if (outlin=="LogR") {
        X <- as.matrix(Idt@LogR)
        Vind <- (p+1):(2*p)
      }
    }
    rownames(X) <- Idt@ObsNames
    k <- robustbase::h.alpha.n(alpha,n,p)
    if (CorrF=="smallsmp") dhn <- robustbase::.MCDcons(p,k/n)*robustbase::.MCDcnp2(p,n,alpha)
    else if (CorrF=="consistent") dhn <- robustbase::.MCDcons(p,k/n)
    else if (CorrF=="none") dhn <- 1.

    if (!force) {
      maxnCk <- 10000000
      nCk <- choose(p,k) 
      if (nCk> maxnCk) stop(paste("fulltle might take too long since",nCk,"different subsets",
                        "need to be evaluated.\nTo proceed anyway set the 'force'",
                        "argument to TRUE, otherwise try the fasttle method instead.\n"))
    }                    

    c0 <- -n*(p*(log(2*pi)+1))/2
    bestCrt <- Inf
    if (SelCrit=="BIC") penC <- log(n)
    else if (SelCrit=="AIC") penC <- 2
    if (outlin=="MidPandLogR")
    {
      for (Cnf in Config) {
        if (Cnf==2) {
          Cftmpsol <- Rfulltle(Idt,k,2,SelCrit,TRUE,...)
        }  else {
          Cftmpsol <- .Call("Cfulltle", X, n, p, k, Cnf, c0, PACKAGE = "MAINT.Data" )
        }
        CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(Cnf,p,p/2)
        Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
        if (CmpCrt < bestCrt) {
          bestCrt <- CmpCrt
          bestSet <- Cftmpsol$Set
        }
      }
    } else {  
      Cftmpsol <- .Call("Cfulltle", X, n, p, k, 1, c0, PACKAGE = "MAINT.Data" )
      CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(1,2*p,p)
      Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
      if (CmpCrt < bestCrt) {
        bestCrt <- CmpCrt
        bestSet <- Cftmpsol$Set
      }
    }
    bestsol <- IdtNmle(Idt[bestSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)

    if (!reweighted && CorrF!="none")
    {
      if (outlin=="MidPandLogR")
      {
        for (Cnf in Config) 
          bestsol@CovConfCases[[Cnf]]$mleSigE <- dhn * bestsol@CovConfCases[[Cnf]]$mleSigE
      } else {
        bestsol@CovConfCases[[1]]$mleSigE <- dhn * bestsol@CovConfCases[[1]]$mleSigE[Vind,Vind]
      }
    }
    if (reweighted) 
    {
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
      bestsol <- IdtNmle(Idt[RewghtdSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
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


    if (otpType=="OnlyEst") 
    { 
      return(finalsol)
    } 
    if (otpType=="SetMD2andEst")
    {    
      return(list(sol=finalsol,rawSet=sort(bestSet),RewghtdSet=RewghtdSet,RobMD2=MD2))
    }
  }
)

