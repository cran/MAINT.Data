setMethod("show",
  signature(object = "IdtOutl"),
  function(object) {
    print(object@outliers)
    invisible(object)
  }
)  

setMethod("plot",
  signature(x = "IdtOutl",y = "missing"),
  function(x,scale=c("linear","log10","log"),RefDist=getRefDist(x),eta=geteta(x),multiCmpCor=getmultiCmpCor(x),...)
  {
    require(ggplot2)   # To be completed (ex: add F and Bets cut-offs and tested ...

    scale <- match.arg(scale)
    MD2 <- getMahaD2(x)
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
      plt <- ggplot(NULL,aes(x=names(MD2),y=MD2)) + geom_point() + geom_abline(intercept=delta,slope=0.)
    } else {
      RewSet <- x@boolRewind
      plt <- ggplot(NULL,aes(x=names(MD2),y=MD2,colour=RewSet)) + geom_point() + 
        geom_abline(intercept=delta1,slope=0.,colour="blue") + geom_abline(intercept=delta2,slope=0.,colour="red")
    }
    if (scale=="linear") {
      plt <- plt + ggtitle("Robust Mahalanobis Distances") + theme(plot.title=element_text(hjust=0.5)) + theme(axis.text.x=element_text(angle=90,hjust=1))       
    } else if (scale=="log10") {
      plt <- plt + ggtitle("Robust Mahalanobis Distances (log10 scale)") +  theme(plot.title=element_text(hjust=0.5)) + 
        theme(axis.text.x=element_text(angle=90,hjust=1)) + coord_trans(y="log10")  
    } else if (scale=="log") {
      plt <- plt + ggtitle("Robust Mahalanobis Distances (logarithmic scale)") + theme(plot.title=element_text(hjust=0.5)) + 
        theme(axis.text.x=element_text(angle=90,hjust=1)) + coord_trans(y="log") 
    }

    plt
  }  
) 

setMethod("getMahaD2",signature(IdtOtl = "IdtOutl"),function(IdtOtl) IdtOtl@MD2) 
setMethod("geteta", signature(IdtOtl = "IdtOutl"), function(IdtOtl) IdtOtl@eta)
setMethod("getRefDist", signature(IdtOtl ="IdtOutl"), function(IdtOtl) IdtOtl@RefDist)
setMethod("getmultiCmpCor", signature(IdtOtl ="IdtOutl"), function(IdtOtl) IdtOtl@multiCmpCor)

