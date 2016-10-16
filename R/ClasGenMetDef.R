setClass("IData",slots=c(MidP="data.frame",LogR="data.frame",ObsNames="character",VarNames="character",
  NObs="numeric",NIVar="numeric"))
setClass("IdtE",slots=c(ModelNames="character",ModelType="character",ModelConfig="numeric",NIVar="numeric",SelCrit="character",
  logLiks="numeric",BICs="numeric",AICs="numeric",BestModel="numeric",SngD="logical"),contains="VIRTUAL")
setClass("IdtSngDE",contains=c("IdtE","VIRTUAL"))
setClass("IdtMxE",slots=c(grouping="factor",Ngrps="numeric"),contains=c("IdtE","VIRTUAL"))
setClass("IdtSngNDE",slots=c(mleNmuE="numeric",mleNmuEse="numeric",CovConfCases="list"),contains="IdtSngDE")
setClass("IdtSngNDRE",slots=c(RobNmuE="numeric",CovConfCases="list"),contains="IdtSngDE")
setClassUnion("extmatrix",c("matrix","NULL"))
setClass("IdtMxNDE",slots=c(Hmcdt="logical",mleNmuE="matrix",mleNmuEse="extmatrix",CovConfCases="list"),contains="IdtMxE")
setClass("IdtMxNDRE",slots=c(Hmcdt="logical",RobNmuE="matrix",CovConfCases="list"),contains="IdtMxE")
setClassUnion("IdtMxtNDE",c("IdtMxNDE","IdtMxNDRE"))
setClassUnion("IdtNDE",c("IdtSngNDE","IdtSngNDRE","IdtMxNDE","IdtMxNDRE"))
setClass("LRTest",slots=c(QuiSq="numeric",df="numeric",pvalue="numeric",H0logLik="numeric",H1logLik="numeric"))
setClass("ConfTests",slots=c(TestRes="list",RestModels="character",FullModels="character"))
setClass("IdtMANOVA",slots=c(NIVar="numeric",grouping="factor",H0res="IdtSngDE",H1res="IdtMxE"),contains="LRTest")
setClass("IdtClMANOVA",contains="IdtMANOVA")
setClass("IdtHetNMANOVA",contains="IdtMANOVA")
setClass("IdtLocSNMANOVA",contains="IdtMANOVA")
setClass("IdtGenSNMANOVA",contains="IdtMANOVA")
setClass("IdtLocNSNMANOVA",contains="IdtMANOVA")
setClass("IdtGenNSNMANOVA",contains="IdtMANOVA")
setClass("Idtlda",slots=c(prior="numeric",means="matrix",scaling="matrix",N="numeric",CovCase="numeric"))
setClass("Idtqda",slots=c(prior="numeric",means="matrix",scaling="array",ldet="numeric",lev="character",CovCase="numeric"))
setClass("IdtSNlocda",slots=c(prior="numeric",ksi="matrix",eta="numeric",scaling="matrix",
  mu="matrix",gamma1="numeric",N="numeric",CovCase="numeric"))
setClass("IdtSNgenda",slots=c(prior="numeric",ksi="matrix",eta="matrix",scaling="array",ldet="numeric",lev="character",
  mu="matrix",gamma1="matrix",CovCase="numeric"))
setClass("IdtSngSNDE",slots=c(CovConfCases="list"),contains="IdtSngDE")
setClass("IdtMxSNDE",slots=c(Hmcdt="logical",CovConfCases="list"),contains="IdtMxE")
setClassUnion("IdtSNDE",c("IdtSngSNDE","IdtMxSNDE"))
setClass("IdtSngNandSNDE",slots=c(NMod="IdtSngNDE",SNMod="IdtSngSNDE"),contains="IdtSngDE")
setClass("IdtMxNandSNDE",slots=c(NMod="IdtMxNDE",SNMod="IdtMxSNDE"),contains="IdtMxE")
setClassUnion("IdtNandSNDE",c("IdtSngNandSNDE","IdtMxNandSNDE"))

setClass("RobEstControl",
  slots=c(
    ncsteps="numeric",
    getalpha="character",
		getkdblstar="character",
    outlin="character",
    trialmethod="character",
    m="numeric",
    reweighted="logical",
    otpType="character"
  ),
  prototype = list(
    ncsteps=200,
    getalpha = "TwoStep75",
    getkdblstar="Twopplusone",
    outlin="MidPandLogR",
    trialmethod="simple",
    m=1,
    reweighted=TRUE,
    otpType="OnlyEst"
  ),
  contains="CovControlMcd"
)

setGeneric("nrow")
setGeneric("ncol")
setGeneric("summary",signature="object")
setGeneric("head",package="utils",signature="x")
setGeneric("tail",package="utils",signature="x")
setGeneric("coef",package="stats",signature="object")
setGeneric("stdEr",package="miscTools",signature="x")
setGeneric("vcov",package="stats",signature="object")
setGeneric("predict",package="stats",signature="object")
setGeneric("lda",package="MASS",signature="x")
setGeneric("qda",package="MASS",signature="x")
setGeneric("mle",
  function(Idt, Model="Normal", CovCase=1:4, SelCrit=c("BIC","AIC"), OptCntrl=list() ,... )
  standardGeneric("mle"))
setGeneric("MANOVA",function(Idt, grouping, Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4, SelCrit=c("BIC","AIC"), 
      Mxt=c("Hom","Het","Loc","Gen"), CVtol=1.0e-5, OptCntrl=list(), onerror=c("stop","warning","silentNull"), ...)
  standardGeneric("MANOVA"))
setGeneric("BestModel",function(ModE,SelCrit=c("IdtCrt","BIC","AIC"))  standardGeneric("BestModel"))
setGeneric("CovCase",function(object)  standardGeneric("CovCase"))
setGeneric("testMod", function(ModE, RestMod=ModE@ModelConfig[2]:length(ModE@ModelConfig), FullMod="Next")
  standardGeneric("testMod"))
setGeneric("H1res",  function(object) standardGeneric("H1res"))
setGeneric("H0res",  function(object) standardGeneric("H0res"))
setGeneric("snda",function(x, grouping, prior="proportions", ...) standardGeneric("snda"))
setGeneric("ObsLogLiks",function(object,Idt,Conf=object@BestModel) standardGeneric("ObsLogLiks"))

setGeneric("fulltle",
  function(Idt, alpha=0.75, reweighted=TRUE, CorrF=c("smallsmp","consistent","none"),
    outlin=c("MidPandLogR","MidP","LogR"),
    CovCase=1:4, SelCrit=c("BIC","AIC"), force=FALSE, otpType=c("OnlyEst","SetMD2andEst"), ... )
  standardGeneric("fulltle"))
setGeneric("fasttle",
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
  standardGeneric("fasttle"))
setGeneric("RobMxtDEst",
  function(Idt,grouping,Mxt=c("Hom","Het"),CovEstMet=c("Pooled","Globdev"),
    CovCase=1:4,SelCrit=c("BIC","AIC"),Robcontrol=RobEstControl(), l1medpar=NULL, ...)
  standardGeneric("RobMxtDEst"))
setGeneric("Roblda",
  function(x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE,  CovEstMet=c("Pooled","Globdev"), SngDMet=c("fasttle","fulltle"),
    Robcontrol=RobEstControl(), ...)
  standardGeneric("Roblda"))
setGeneric("Robqda",
  function(x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, SngDMet=c("fasttle","fulltle"),
      Robcontrol=RobEstControl(), ...) 
  standardGeneric("Robqda"))