setMethod("show",
  signature(object = "IdtOutl"),
  function(object) {
    print(object@outliers)
    cat("MD2\n");
    print(object@MD2[object@outliers])
    invisible(object)
  }
)  

setMethod("plot",
  signature(x = "IdtOutl",y = "missing"),
  function(x,scale=c("linear","log"),RefDist=getRefDist(x),eta=geteta(x),multiCmpCor=getmultiCmpCor(x),
           Obsnames=TRUE,sqrddist=TRUE, ...)
  {
    scale <- match.arg(scale)
    MD2 <- getMahaD2(x)
    if (sqrddist) { y <- MD2 ; ylb <- "MD2" ; trshld <- x@trshld }
    else { y <- sqrt(MD2) ; ylb <- "MD" ; trshld <- sqrt(x@trshld) }

    n <- x@NObs
    p <- x@p
    h <- x@h

    oneminuseta <- 1-eta
    if (multiCmpCor=="always") {
      oneminusalpha <- oneminuseta^(1/n)
    } else if (multiCmpCor=="never") {
      oneminusalpha <- oneminuseta
    }

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

   if (RefDist!="CerioliBetaF")  {
     if (scale=="linear") {
       plot.default(y,main="Robust Mahalanobis Distances",xlab="",xaxt="n",ylab=ylb,...)
     } else if (scale=="log") {
       plot.default(y,main="Robust Mahalanobis Distances (log scale)",xlab="",xaxt="n",ylab=paste("ln",ylb),log="y",...)
     }
    if (Obsnames) axis(1,1:n,labels=names(x@MD2),las=2,...)
    abline(h=trshld[1])
    if (length(trshld)==2) abline(h=trshld[2]) 
   } else {
      RewSet <- x@boolRewind
      if (scale=="linear") {
        plot.default(y,main="Robust Mahalanobis Distances",xlab="",xaxt="n",type="n",ylab=ylb,...)
      } else if (scale=="log") {
        plot.default(MD2,main="Robust Mahalanobis Distances (log scale)",xlab="",xaxt="n",type="n",log="y",ylab=paste("ln",ylb),...)
      }
      if (Obsnames) axis(1,1:n,labels=names(x@MD2),las=2,...)
      abline(h=trshld[1])
      if (length(trshld)==2) abline(h=trshld[2]) 
      RewSetInd <- which(RewSet==TRUE)
      UnRewSetInd <- which(RewSet==FALSE)
      points(x=RewSetInd,y=MD2[RewSet],pch=19,col="blue")
      points(x=UnRewSetInd,y=MD2[!RewSet],pch=19,col="red")
      abline(h=delta1,col="blue")
      abline(h=delta2,col="red")
   }
  }
)


setMethod("getMahaD2",signature(SdtOtl = "IdtOutl"),function(SdtOtl) SdtOtl@MD2) 
setMethod("geteta", signature(SdtOtl = "IdtOutl"), function(SdtOtl) SdtOtl@eta)
setMethod("getRefDist", signature(SdtOtl ="IdtOutl"), function(SdtOtl) SdtOtl@RefDist)
setMethod("getmultiCmpCor", signature(SdtOtl ="IdtOutl"), function(SdtOtl) SdtOtl@multiCmpCor)

