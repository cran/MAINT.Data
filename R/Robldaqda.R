setMethod("Roblda",
  signature(x = "IData"),
  function(x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE,  CovEstMet=c("Pooled","Globdev"), SngDMet=c("fasttle","fulltle"),
    Robcontrol=RobEstControl(), ...)
  {
    CovEstMet <- match.arg(CovEstMet)
    SelCrit <- match.arg(SelCrit)
    SngDMet <- match.arg(SngDMet)
    if (SngDMet=="fulltle")  { stop("Roblda with fulltle not implemented yet.\n") }

    if (length(subset) < x@NObs)
    {
      x <- x[subset,]
      grouping <- grouping[subset]
    }
    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    n <- x@NObs
    p <- 2*x@NIVar
    k <- length(grplvls)
    if (length(grplvls)==1)
    { 
      errmsg <- "The data belongs to one single group. A partition into at least two different groups is required\n"
      if (silent==FALSE) {
        stop(errmsg)
      }  else {
        warning(errmsg)
        return(NULL)
      }
    }
    MxtDEst <- RobMxtDEst(x,grouping,Mxt="Hom",CovCase=CovCase,SelCrit=SelCrit,CovEstMet=CovEstMet,...)
    nk <- as.numeric(table(grouping))
    if (sum(nk)!=n)  { stop("Dimensions of the x and grouping arguments do not agree with each other\n") }
    glbmeans <- colMeans(cbind(x@MidP,x@LogR))
    mugdev <- scale(MxtDEst@RobNmuE,center=glbmeans,scale=FALSE)
    vnames <- unlist(dimnames(MxtDEst@RobNmuE)[2]) 
    B <- matrix(0.,nrow=p,ncol=p,dimnames=list(vnames,vnames))
    for (g in 1:k) B <- B + (nk[g]/n) * outer(mugdev[g,],mugdev[g,]) 
    selmodel <- BestModel(MxtDEst)

    Ilda(Conf=selmodel,p=p,nk=nk,prior=prior,means=MxtDEst@RobNmuE,W=MxtDEst@CovConfCases[[selmodel]]$RobSigE,B=B,egvtol=egvtol,...)
  }
)

setMethod("Robqda",
  signature(x = "IData"),
  function(x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, SngDMet=c("fasttle","fulltle"),
      Robcontrol=RobEstControl(), ...) 
  {
    SelCrit <- match.arg(SelCrit)
    SngDMet <- match.arg(SngDMet)

    if (length(subset) < x@NObs)
    {
      x <- x[subset,]
      grouping <- grouping[subset]
    }
    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    if (length(grplvls)==1)
    { 
      errmsg <- "The data belongs to one single group. A partition into at least two different groups is required\n"
      if (silent==FALSE) {
        stop(errmsg)
      }  else {
        warning(errmsg)
        return(NULL)
      }
    }
    MxtDEst <- RobMxtDEst(x,grouping,Mxt="Het",CovCase=CovCase,SelCrit=SelCrit,...)
    selmodel <- BestModel(MxtDEst)

    Iqda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(grouping)),lev=grplvls,
      prior=prior,means=MxtDEst@RobNmuE,Wg=coef(MxtDEst,selmodel)$Sigma)
  }
)

