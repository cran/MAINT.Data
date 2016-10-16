setMethod("summary",
  signature(object = "IData"),
  function (object) 
  {
    Rngsumar <- summary(exp(object@LogR))
    dimnames(Rngsumar)[[2]] <- paste(object@VarNames,".Range",sep="")
    res <- list(MidPsumar=summary(object@MidP),Rngsumar=Rngsumar,LogRsumar=summary(object@LogR))
    class(res) <- "summaryIData"
    res
  }
)

setMethod("show",
  signature(object = "IData"),
  function (object) 
  {
    printrow <- function(Bnds,NIVar) 
    { 
      cat(Bnds[1],"  ")
      for (j in 2:(NIVar+1))
        cat("[",format(Bnds[j],width=8,digits=5,justify="centre"),", ",
          format(Bnds[NIVar+j],width=8,digits=5,justify="centre"),"]  ",sep="") 
      cat("\n")
    }

    HalfRange <- exp(object@LogR)/2
    LB <- object@MidP - HalfRange
    UB <- object@MidP + HalfRange
    lobsname <- max(nchar(object@ObsNames))
    centerp <- lstpos <- 0
    cat <- "      "
    flength <- nchar(format(LB[1,1],width=8,digits=5)) + nchar(format(UB[1,1],width=8,digits=5)) + 7 
    for (j in 1:object@NIVar) {
      centerp <- centerp + flength
      hlfname <- nchar(object@VarNames[j])/2
      cat(rep(" ",centerp-lstpos-hlfname-1),object@VarNames[j],sep="" )
      lstpos <- centerp + hlfname
    }
    cat("\n") 
    apply(cbind(format(object@ObsNames,width=lobsname),LB,UB),1,printrow,NIVar=object@NIVar)
    invisible(object)
  }
)

setMethod("nrow",signature(x = "IData"),function (x) {x@NObs})

setMethod("ncol",signature(x = "IData"),function (x) {x@NIVar})

setMethod("rownames",signature(x = "IData"),function (x) {x@ObsNames})

setMethod("colnames",signature(x = "IData"),function (x) {x@VarNames})


setMethod("head",
  signature(x = "IData"),
  function (x,n=6L) 
  {
    if (n>0) {
      x[1:n,] 
    } else {
      x[-1:n,]
    }
  }
)

setMethod("tail",
  signature(x = "IData"),
  function (x,n=6L) 
  {
    if (n>0) {
      x[(x@NObs-n+1):x@NObs,] 
    } else {
      x[(-x@NObs-n-1):-x@NObs,]
    }
  }
)

print.summaryIData <- function(x,...) 
{
  cat("Mid-Points summary:\n") ; print(x$MidPsumar)
  cat("Ranges summary:\n") ; print(x$Rngsumar)
  cat("Log-Ranges summary:\n") ; print(x$LogRsumar)
  invisible(x)
}

IData <- function(Data,
  Seq=c("LbUb_VarbyVar","MidPLogR_VarbyVar","AllLb_AllUb","AllMidP_AllLogR"),VarNames=NULL,ObsNames=row.names(Data))
{
  if (!(is.data.frame(Data))) stop("First argument of IData must be a data frame\n")
  p <- ncol(Data)  # Total number of Interval variable bounds
  q <- p/2	 # Number of Interval variables
  if (floor(q) != q) stop("Number of columns of Data ( =",p,") must be an even number\n")
  Seq <- match.arg(Seq)
  if (  (Seq == "LbUb_VarbyVar") || (Seq == "AllLb_AllUb") )
  {
    if (Seq == "LbUb_VarbyVar") { Lbnd <- Data[,2*(0:(q-1))+1] ; Ubnd <- Data[,2*(1:q)] }
    if (Seq == "AllLb_AllUb")   { Lbnd <- Data[,1:q] ; Ubnd <- Data[,(q+1):p] }
    MidP <- (Lbnd+Ubnd)/2
    LogR <- log(Ubnd-Lbnd)
    if (any(is.na(LogR))) stop("Invalid data")
  } else {
    if (Seq == "MidPLogR_VarbyVar") { MidP <- Data[,2*(0:(q-1))+1] ; LogR <- Data[,2*(1:q)] }
    if (Seq == "AllMidP_AllLogR")   { MidP <- Data[,1:q] ; LogR <- Data[,(q+1):p] }
  }
  if (is.null(VarNames)) VarNames <- paste("I",1:q,sep="")
  if (!is.data.frame(MidP)) MidP <- as.data.frame(MidP)
  if (!is.data.frame(LogR)) LogR <- as.data.frame(LogR)
  names(MidP) <- paste(VarNames,".MidP",sep="")
  names(LogR) <- paste(VarNames,".LogR",sep="")

  new("IData",MidP=MidP,LogR=LogR,ObsNames=ObsNames,VarNames=VarNames,NObs=nrow(MidP),NIVar=q)
}

# Standard operators for IData objects

# Indexing and assignement

"[.IData" <- function(x,rowi=1:x@NObs,coli=1:x@NIVar,...)
{
  IData(cbind(x@MidP[rowi,coli,drop=FALSE],x@LogR[rowi,coli,drop=FALSE]),
    Seq="AllMidP_AllLogR",VarNames=x@VarNames[coli],ObsNames=x@ObsNames[rowi])
}

"[<-.IData" <- function(x,rowi=1:x@NObs,coli=1:x@NIVar,value)
{
  if (!is(value,"IData")) stop("Argument value is not an IData object\n")
  x@MidP[rowi,coli] <- value@MidP
  x@LogR[rowi,coli] <- value@LogR
  x
}

# Comparison

"==.IData" <- function(x,y)
{
  CompIvalue <- function(Ival) Ival[1,1] == Ival[1,2] && Ival[2,1] == Ival[2,2]

  if (!is(y,"IData")) 
    stop("Trying to compare an IData object with an object of a diferent type\n")
  if ( x@NObs != y@NObs || x@NIVar != y@NIVar )
    stop("== only defined for equally-sized IData objects\n")
  TmpArray <- array(dim=c(x@NObs,x@NIVar,2,2))
  for (j in 1:x@NIVar)  {
    TmpArray[,j,,1] <- cbind(x@MidP[,j],x@LogR[,j])
    TmpArray[,j,,2] <- cbind(y@MidP[,j],y@LogR[,j])
  }
  apply(TmpArray,c(1,2),CompIvalue)
}

"!=.IData" <- function(x,y)  
{
  if (!is(y,"IData"))
    stop("Trying to compare an IData object with an object of a diferent type\n")
  if ( x@NObs != y@NObs || x@NIVar != y@NIVar )
    stop("!= only defined for equally-sized IData objects\n")
  !(x==y)
}


