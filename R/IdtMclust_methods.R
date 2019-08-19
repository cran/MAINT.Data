setMethod("pcoordplot",
  signature(x = "IdtMclust"),
  function(x,title="Parallel Coordinate Plot",Seq=c("AllMidP_AllLogR","MidPLogR_VarbyVar"),G=BestG(x),...)
  {
    if (requireNamespace("GGally",quietly=TRUE)==FALSE) {
      stop("Required package GGally is not installed\n")
    }

    q <- x@NIVar   # Numbef or Interval-valued variables
    p <- 2*q       # Total number of MidPoints and LogRanges
    
    Seq <- match.arg(Seq)
    if (Seq == "AllMidP_AllLogR")  {
      DF <- as.data.frame(t(mean(x)))
    }  else if (Seq == "MidPLogR_VarbyVar")  {
      DF <- as.data.frame(t(mean(x)[rep(1:q,each=2)+rep(c(0,q),q),]))
    }
    DF <- cbind(DF,Component=paste("CP",1:G,sep=""))
    plt <- ggparcoord(DF, columns = 1:p, groupColumn = 'Component') + ggtitle(title) +
      theme_minimal() + theme(plot.title=element_text(hjust=0.5)) + theme(axis.text.x=element_text(angle=90,hjust=1)) 
       
    plt <- plt + theme(axis.title = element_blank())

    plt <- plt + theme(axis.text.y = element_blank())
    plt <- plt + theme(panel.grid.major.y = element_blank())

    plt <- plt + theme(panel.grid.minor = element_blank())
    plt <- plt + theme(panel.grid.major.x = element_line(color = "#bbbbbb"))
    min_y <- min(plt$data$value)
    max_y <- max(plt$data$value)
    pad_y <- (max_y - min_y) * 0.1
    lab_x <- rep(1:p, times = 2) # 2 times, 1 for min 1 for max
    lab_y <- rep(c(min_y - pad_y, max_y + pad_y), each = p)

    z <- c(sapply(DF[,1:p], min), sapply(DF[, 1:p], max))
    absz <- abs(z)
    z_scinot <- which(absz<1e-3|absz>=1e4)
    z_3dig <- which(absz>=1e-3&absz<1.)
    z_2dig <- which(absz>=1.&absz<10.)
    z_1dig <- which(absz>=10.&absz<100.)
    z_0dig <- which(absz>=100.&absz<1e4)
    
    lab_z <- character(length(z))
    if (length(z_scinot)>0) {
      lab_z[z_scinot] <- formatC(z[z_scinot],digits=1,format="e")
    }
    if (length(z_3dig)>0) {
      lab_z[z_3dig] <- formatC(z[z_3dig],digits=3,format="f")
    }
    if (length(z_2dig)>0) {
      lab_z[z_2dig] <- formatC(z[z_2dig],digits=2,format="f")
    }
    if (length(z_1dig)>0) {
      lab_z[z_1dig] <- formatC(z[z_1dig],digits=1,format="f")
    }
    if (length(z_0dig)>0) {
      lab_z[z_0dig] <- formatC(z[z_0dig],digits=0,format="f")
    }

    plt <- plt + annotate("text", x = lab_x, y = lab_y, label = lab_z, size = 3)
    
    print(plt)
  }
)

setMethod("show",
  signature(object = "IdtMclust"),
  function(object)
  {
    cat("\'", class(object)[1], "\' model object:\n", sep = "")
    if (object@Hmcdt) {
      str1 <- "Homoscedastic"
    } else { 
      str1 <- "Heteroscedastic"
    }
    cat(" best model: ", str1, " setup with covariance configuration C",object@BestC," and ", object@BestG, " components\n",sep="") 
    invisible()
  }
)  

setMethod("summary",
  signature(object = "IdtMclust"),
  function(object, parameters = FALSE, classification = FALSE, ...)
  {
    G  <- object@BestG
    pro <- object@parameters$pro
    if(is.null(pro)) pro <- 1
    names(pro) <- paste("CP",seq(G),sep="")
    mean <- object@parameters$mean
    covariance <- object@parameters$covariance
    title <- paste("Gaussian finite mixture model fitted by EM algorithm")
    if (object@Hmcdt) {
      str1 <- "Homecedastic"
    } else { 
      str1 <- "Heteroscedastic"
    }
    modelName <- paste(str1, " C",object@BestC,sep="") 
    obj <- list(
      title = title, modelName =modelName,  Hmcdt = object@Hmcdt,     
      NObs = object@NObs, NIVar = object@NIVar, G = G,  
      loglik = object@logLik, bic = object@bic,
      pro = pro, mean = mean, covariance = covariance,
      classification = object@classification,
      printParameters = parameters, printClassification = classification
    )
    class(obj) <- "summaryIdtMclust"
    return(obj)
  }
) 

print.summaryIdtMclust <- function(x, ...)
{
   n <- x$NObs
   p <- 2*x$NIVar
   cat(rep("-", nchar(x$title)),"\n",sep="")
   cat(x$title, "\n")
   cat(rep("-", nchar(x$title)),"\n",sep="")
   cat(x$modelName," model with ", x$G, ifelse(x$G > 1, " components\n", " component\n"),sep = "") 
   tab <- data.frame("log-likelihood" = x$loglik, "NObs" = n, "BIC" = x$bic, row.names = "")
   print(tab)
   cat("\nClustering table:")
   if (x$G==1) {
     cat("\nCP1\n",n,"\n")
   } else {
     print(table(factor(x$classification, levels = paste("CP",seq_len(x$G),sep=""))))
   }
   if(x$printParameters) {
     cat("\nMixing probabilities:\n")
     print(x$pro)
     cat("\nMeans:\n")
     print(x$mean)
     cat("\nStandard deviations:\n")
     if (x$Hmcdt) {
       print(sqrt(diag(x$covariance[,,1])))
     } else {
       stdv <- matrix(nrow=p,ncol=x$G,dimnames=list(rownames(x$mean),colnames(x$mean)))
       for(g in 1:x$G) stdv[,g] <- sqrt(diag(x$covariance[,,g]))
       print(stdv)
     }
     cat("\nCorrelations:\n")
     if (x$Hmcdt) {
       print(cov2cor(x$covariance[,,1]))
     } else { 
       for(g in 1:x$G) { 
         cat("[,,CP", g, "]\n", sep = "")
         print(cov2cor(x$covariance[,,g])) 
       }
     }
   }
   if(x$printClassification) {
     cat("\nClassification:\n")
     print(x$classification)
   }
   invisible(x)
 }

#Accessor methods

setMethod("parameters",signature(x = "IdtMclust"),function(x) x@parameters)
setMethod("pro",signature(x = "IdtMclust"),function(x) x@parameters$pro)
setMethod("mean",signature(x = "IdtMclust"),function(x) x@parameters$mean)
setMethod("var",signature(x = "IdtMclust"),function(x) x@parameters$covariance)
setMethod("classification",signature(x = "IdtMclust"),function(x) x@classification)

setMethod("SelCrit",signature(x = "IdtMclust"),function (x) x@SelCrit)
setMethod("Hmcdt",signature(x = "IdtMclust"),function (x) x@Hmcdt)
setMethod("BestG",signature(x = "IdtMclust"),function (x) x@BestG)
setMethod("BestC",signature(x = "IdtMclust"),function (x) x@BestC)
setMethod("PostProb",signature(x = "IdtMclust"),function(x) x@z)

setMethod("logLik",signature(x = "IdtMclust"),function(x) x@logLik)
setMethod("BIC",signature(x = "IdtMclust"),function(x) x@bic)
setMethod("AIC",signature(x = "IdtMclust"),function(x) x@aic)

setMethod("cor",signature(x ="IdtMclust"),
  function(x)
 {
    if (x@Hmcdt==TRUE) {
      ncov <- 1 
    } else {
      ncov <- x@BestG
    }

    covdim <- 2*x@NIVar
    covnames <- dimnames(x@parameters$covariance)[[1]]   
    res <- array(dim=c(covdim,covdim,ncov),dimnames=list(covnames,covnames,NULL))
    for (g in 1:ncov) res[,,g] <- cov2cor(x@parameters$covariance[,,g])

    res
  }
) 



