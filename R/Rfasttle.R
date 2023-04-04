trialstep <- function(Sdt,n,p,c0,kdblstar,maxrefstps,method=c("simple","Poolm"),m=NULL,
  Config,SelCrit,convrg=.Machine$double.eps,OptCntrl=list())
{
  method <- match.arg(method)
  if (method=="simple")  {
    IdtNmle(Sdt[sample.int(n,kdblstar),],OptCntrl=OptCntrl,CovCaseArg=FALSE,Config=Config,SelCrit=SelCrit)
  }  else if (method=="Poolm") {
    if (is.null(m))  { stop("m argument for trial step was not specified\n") }
    indlim <- floor((0:m)*n/m)
    permut <- sort.int(runif(n),index.return=TRUE)$ix
    Set <- NULL
    for (i in 1:m) {
      sbstind <- permut[(indlim[i]+1):indlim[i+1]]
      sbstIdt <- Sdt[sbstind,]
      initsmpl <- sample.int(sbstIdt@NObs,kdblstar)
      tmpsol <- IdtNmle(sbstIdt[initsmpl,],OptCntrl=OptCntrl,CovCaseArg=FALSE,Config=Config,SelCrit=SelCrit)
      kstar <- ceiling((sbstIdt@NObs+p+1)/2)
      tmpsol <- refinementstep(tmpsol,sbstIdt,p,c0,kstar,maxrefstps,Config,SelCrit,convrg,OptCntrl)
      Set <- c(Set,tmpsol$Set)
    }
    IdtNmle(Sdt[Set,],OptCntrl=OptCntrl,CovCaseArg=FALSE,Config=Config,SelCrit=SelCrit)
  }
}

refinementstep <- function(tmpsol,Sdt,p,c0,k,maxnsteps,Config,SelCrit,convrg=.Machine$double.eps,OptCntrl=list())
{
  Cftmpsol <- vector("list",5)
  CfSet <- matrix(nrow=5,ncol=k)
  CFCrit <- rep(NA,5)
  X <- cbind(Sdt@MidP,Sdt@LogR)
  for (Cf in Config)  {
    Cftmpsol[[Cf]] <- tmpsol
    for (step in 1:maxnsteps) {
      if (step==1)  {
        prevLogLik <- -Inf 
      }  else  {
        prevLogLik <- Cftmpsol[[Cf]]@logLiks[Cf]
      }
      Xdev <- scale(X,center=Cftmpsol[[Cf]]@mleNmuE,scale=FALSE)
      if (Cf!=5)
      {
        SigISr <- t(backsolve(chol(Cftmpsol[[Cf]]@CovConfCases[[Cf]]$mleSigE),diag(p)))
        c1 <- sum(log(diag(SigISr)))
        obsLogL <- apply(Xdev,1,ILogLikNC1,SigmaSrInv=SigISr,const=c0+c1)
      }  else  {
        IVar <- 1./diag(Cftmpsol[[Cf]]@CovConfCases[[5]]$mleSigE)
        c1 <- -0.5*sum(log(IVar))
        obsLogL <- apply(Xdev,1,ILogLikDNC,IVar=IVar,const=c0+c1)
      }
      CfSet[Cf,] <- sort(obsLogL,decreasing=TRUE,index.return=TRUE)$ix[1:k]
      Cftmpsol[[Cf]] <- IdtNmle(Sdt[CfSet[Cf,],],OptCntrl=OptCntrl,CovCaseArg=FALSE,Config=Cf,SelCrit=SelCrit)
      if (  is.finite(prevLogLik) && (Cftmpsol[[Cf]]@logLiks[Cf]-prevLogLik)/abs(prevLogLik) < convrg  )  {
        break
      }
    }
    if (SelCrit=="BIC")  {
      CFCrit[Cf] <- Cftmpsol[[Cf]]@BICs[Cf] 
    }  else if (SelCrit=="AIC")  {
      CFCrit[Cf] <- Cftmpsol[[Cf]]@AICs[Cf]
    }
  }
  BestCf <- which.min(CFCrit)		
  list(LogLik=Cftmpsol[[BestCf]]@logLiks[BestCf],Set=CfSet[BestCf,])
}

Rfasttle <- function(Sdt,kdblstar=2*Sdt@NIVar+1,k=ceiling((Sdt@NObs+2*Sdt@NIVar+1)/2),nrep=500,
  Config=2,SelCrit=c("BIC","AIC"),maxrefstps=100,trialmethod=c("simple","Poolm"),m=NULL,...)
{
  SelCrit <-  match.arg(SelCrit)
  trialmethod <- match.arg(trialmethod)
  n <- Sdt@NObs	
  p <- 2*Sdt@NIVar
  c0 <- -0.5*(p*log(2*pi))
  bestsol <- NULL
  for (rep in 1:nrep)  {
    trialsol <- trialstep(Sdt,n,p,c0,kdblstar,maxrefstps,method=trialmethod,m=m,Config=Config,SelCrit=SelCrit,...)
    tmpsol <- refinementstep(trialsol,Sdt,p,c0,k,maxrefstps,Config,SelCrit,...)
    if (is.null(bestsol) || tmpsol$LogLik > bestsol$LogLik)  {
      bestsol <- tmpsol
    }
  }

  finalsol <- IdtNmle(Sdt[bestsol$Set,],CovCaseArg=FALSE,Config=Config,SelCrit=SelCrit,...)
  BestModel <- finalsol@BestModel

  list(LogLik=finalsol@logLiks[BestModel],Set=bestsol$Set,raw.cov=finalsol@CovConfCases[BestModel]$mleSigE)
}
