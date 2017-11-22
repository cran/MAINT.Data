setMethod("show",
  signature(object = "IdtMclust"),
  function(object)
  {
    cat("\'", class(object)[1], "\' model object:\n", sep = "")
    if (object@Hmcdt) {
      str1 <- "Homecedastic"
    } else { 
      str1 <- "Heterocedastic"
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
    names(pro) <- seq(G)
    mean <- object@parameters$mean
    covariance <- object@parameters$covariance
    title <- paste("Gaussian finite mixture model fitted by EM algorithm")
    if (object@Hmcdt) {
      str1 <- "Homecedastic"
    } else { 
      str1 <- "Heterocedastic"
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

print.summaryIdtMclust <- function(x, digits = getOption("digits"), ...)
{
  cat(rep("-", nchar(x$title)),"\n",sep="")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)),"\n",sep="")
  cat(x$modelName," model with ", x$G, ifelse(x$G > 1, " components\n", " component\n"),sep = "") 
  tab <- data.frame("log-likelihood" = x$loglik, "NObs" = x$NObs, "BIC" = x$bic, row.names = "")
  print(tab, digits = digits)
  cat("\nClustering table:")
  print(table(factor(x$classification, levels = seq_len(x$G))),digits = digits)
  if(x$printParameters) {
    cat("\nMixing probabilities:\n")
    print(x$pro, digits = digits)
    cat("\nMeans:\n")
    print(x$mean, digits=digits)
    cat("\nStandard deviations:\n")
    if (x$Hmcdt) {
      print(sqrt(diag(x$covariance[,,1])),digits=digits)
    } else {
      stdv <- matrix(nrow=2*x$NIVar,ncol=x$G,dimnames=list(rownames(x$mean),NULL))
      for(g in 1:x$G) stdv[,g] <- sqrt(diag(x$covariance[,,g]))
      print(stdv, digits=digits)
    }
    cat("\nCorrelations:\n")
    if (x$Hmcdt) {
      print(cov2cor(x$covariance[,,1]),digits=digits)
    } else { 
      for(g in 1:x$G) { 
        cat("[,,", g, "]\n", sep = "")
        print(cov2cor(x$covariance[,,g]),digits=digits) 
      }
    }
  }
  if(x$printClassification) {
    cat("\nClassification:\n")
    print(x$classification, digits = digits)
  }
  invisible(x)
}


