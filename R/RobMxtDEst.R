setMethod("RobMxtDEst",
  signature(Idt = "IData"),
  function(Idt,grouping,Mxt=c("Hom","Het"),CovEstMet=c("Pooled","Globdev"),
    CovCase=1:4,SelCrit=c("BIC","AIC"),Robcontrol=RobEstControl(), l1medpar=NULL, ...)
  {
    if (!requireNamespace("robustbase",quietly=TRUE)) 
      stop("fasttle needs the robustbase package to work. Please install it\n")

    Mxt <- match.arg(Mxt)
    CovEstMet <- match.arg(CovEstMet)
    SelCrit <- match.arg(SelCrit)
    n <- Idt@NObs  
    p <- 2*Idt@NIVar  
    q <- p/2
    if (length(grouping)!=n)  {
      stop("The size of the grouping factor does not agree with the number of observations in the data set supplied")
    }

    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    ng <- as.integer(table(grouping))
    ngrps <- length(ng)
    grpRobE <- vector("list",ngrps)
    Xnams <- names(cbind(Idt@MidP,Idt@LogR))
    anams <- list(Xnams,Xnams,levels(grouping))
    RobNmuE <- matrix(nrow=ngrps,ncol=p,dimnames=list(grplvls,Xnams))

    nCovC <- length(CovCase)
    CovConfC <- vector("list",nCovC)
    logLiks <- numeric(nCovC) 
    AICs <- numeric(nCovC) 
    BICs <- numeric(nCovC) 
    names(logLiks) <- names(AICs) <- names(BICs) <- names(CovConfC) <- modnames <- paste("NModCovC",CovCase,sep="")

    if (Mxt=="Hom" && CovEstMet=="Globdev")
    {
      X <- cbind(Idt@MidP,Idt@LogR)
      Xdev <- matrix(nrow=n,ncol=p) 
      Xgl1med <- matrix(nrow=ngrps,ncol=p) 
      if (!is.null(l1medpar)) {
        MaxStep <- ifelse(is.null(l1medpar$MaxStep),200,l1medpar$MaxStep)
        ItTol <- ifelse(is.null(l1medpar$ItTol),10^-8,l1medpar$ItTol)
        trace <- ifelse(is.null(l1medpar$trace),0,l1medpar$trace)
      }
      for (g in 1:ngrps)
      {
        gind <- which(grouping==grplvls[g])
        Xg <- X[gind,]
        if (is.null(l1medpar)) {
          Xgl1med[g,] <- l1median(Xg)
        }  else {
          if (is.null(l1medpar$m.init)) {
            m.init <- robustbase::colMedians(as.matrix(Xg))
          } else {
            m.init <- l1medpar$m.init
          }
          Xgl1med[g,] <- l1median(Xg,MaxStep,ItTol,trace,m.init)
        }
        Xdev[gind,] <- scale(Xg,center=Xgl1med[g,],scale=FALSE)
      }
    }
    alpha <- Robcontrol@alpha
    nSteps <- ifelse(Robcontrol@getalpha=="TwoStep",2,1)
    for (Steps in 1:nSteps)
    {
      for (CovC in 1:nCovC)
      {
        if (Mxt=="Hom")
        {
          CovConfC[[CovC]] <- list(
            RobSigE=matrix(0.,nrow=p,ncol=p,dimnames=list(Xnams,Xnams)),logLik=NULL,AIC=NULL,BIC=NULL
          )
          trmdn <- round(alpha*n)
          if (CovEstMet=="Pooled")  {
            X <- NULL
            for (g in 1:ngrps)  {      # Things to do: add a fulltle option !!!
              Idtg <- Idt[grouping==grplvls[g],]
              grpRobE[[g]] <- fasttle(Idtg,CovC,SelCrit,
                 getalpha="NO",otpType="SetMD2andEst",Robcontrol=Robcontrol, ...)
              RobNmuE[g,] <- grpRobE[[g]]$sol@RobNmuE
              CovConfC[[CovC]]$RobSigE <- 
                CovConfC[[CovC]]$RobSigE + (ng[g]/n) * grpRobE[[g]]$sol@CovConfCases[[CovC]]$RobSigE 
            }
            regset <- ifelse(Robcontrol@reweighted,grpRobE[[g]]$RewghtdSet,grpRobE[[g]]$rawSet)
            X <- rbind(X,cbind(Idtg@MidP[regset,],Idtg@LogR[regset,]))
          }
          else if (CovEstMet=="Globdev") {
            if  (Robcontrol@getkdblstar=="Twopplusone") { 
              kdblstar <- 2*Idt@NIVar+1
            }  else {
              if (!is.finite(Robcontrol@getkdblstar)) {
                stop("Wrong value for Robcontrol parameter getkdblstar\n")
              }
              kdblstar <- Robcontrol@getkdblstar 
            }
            RobE <- fasttle1(Xdev,CovC,SelCrit,alpha,Robcontrol@nsamp,Robcontrol@ncsteps,Robcontrol@trace,
              Robcontrol@use.correction,kdblstar,Robcontrol@outlin,Robcontrol@trialmethod,Robcontrol@m,
              Robcontrol@reweighted,"OnlyEst",Idt@VarNames,...)
            for (g in 1:ngrps)  RobNmuE[g,] <- Xgl1med[g,] + RobE@RobNmuE
            CovConfC[[CovC]]$RobSigE <- RobE@CovConfCases[[CovC]]$RobSigE 
          }
          Xdev <- scale(X[grouping==levels(grouping)[1],],center=RobNmuE[1,],scale=FALSE)
          for (g in 2:ngrps)
            Xdev <- rbind(Xdev,scale(X[grouping==levels(grouping)[g],],center=RobNmuE[g,],scale=FALSE))
          logdet <- pdwt.solve(CovConfC[[CovC]]$RobSigE,silent=TRUE,onlylogdet=TRUE)
          if (is.null(logdet))  {
            logLiks[CovC] <- CovConfC[[CovC]]$logLik <- -Inf
          }  else  {
            logLiks[CovC] <- CovConfC[[CovC]]$logLik <- -trmdn*(p*(log(2*pi)+1)+logdet)/2
          }		
        }  else if (Mxt=="Het") {
          CovConfC[[CovC]] <- list( RobSigE=array(dim=c(p,p,ngrps),dimnames=anams),logLik=NULL,AIC=NULL,BIC=NULL )
          for (g in 1:ngrps) {
            grpRobE[[g]] <- fasttle(Idt[grouping==grplvls[g],],CovC,SelCrit,
              getalpha="NO",otpType="OnlyEst",Robcontrol=Robcontrol, ...)
            RobNmuE[g,] <- grpRobE[[g]]@RobNmuE
            CovConfC[[CovC]]$RobSigE[,,g] <- grpRobE[[g]]@CovConfCases[[CovC]]$RobSigE
          } 
          logLiks[CovC] <- grpRobE[[1]]@logLiks[CovC]
          for (g in 2:ngrps) logLiks[CovC] <- logLiks[CovC] + grpRobE[[g]]@logLiks[CovC]
        }
        CovConfC[[CovC]]$logLik <- logLiks[CovC]
        nmodelfreepar <- npar(CovC,p,q,Ngrps=ngrps,Mxt=Mxt)
        CovConfC[[CovC]]$AIC <- AICs[CovC] <- -2*logLiks[CovC] + 2*nmodelfreepar
        trmdn <- round(alpha*n)
        CovConfC[[CovC]]$BIC <- BICs[CovC] <- -2*logLiks[CovC] + log(trmdn)*nmodelfreepar
      }
      if (SelCrit=="AIC")  {
        bestmod = which.min(AICs)
      }  else if (SelCrit=="BIC")  {
        bestmod = which.min(BICs)
      }

      if(Robcontrol@getalpha=="TwoStep" && Steps==1)
      {
        X <- data.frame(cbind(Idt@MidP,Idt@LogR))
        nOtls <- 0.
        for (g in 1:ngrps)
        {
          if (Mxt=="Hom")  {
            nOtls <- nOtls + MDOtlDet(X[grouping==grplvls[g],],RobNmuE[g,],CovConfC[[bestmod]]$RobSigE,0.025,otp="onlycnt") 
          } else if (Mxt=="Het") {
            nOtls <- nOtls + MDOtlDet(X[grouping==grplvls[g],],RobNmuE[g,],CovConfC[[bestmod]]$RobSigE[,,g],0.025,otp="onlycnt") 
          }
        }
        alpha <- 1.-nOtls/n
      }

    }

    new( "IdtMxNDRE", ModelNames=modnames,ModelType=rep("Normal",nCovC),ModelConfig=1:nCovC,
      grouping=grouping,Hmcdt=(Mxt=="Hom"),RobNmuE=RobNmuE,CovConfCases=CovConfC,
      SelCrit=SelCrit,NIVar=q,logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=ngrps)

  }
)

