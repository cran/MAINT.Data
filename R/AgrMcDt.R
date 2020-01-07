AgrMcDt <- function(MicDtDF,agrby,agrcrt="minmax")
{
  if (!(is.data.frame(MicDtDF))) stop("First argument of AgMicroData must be a data frame\n")
#  if (class(agrby)!="factor") stop("Argument agrby is not a factor\n")
  if (class(agrby)[1]!="factor") stop("Argument agrby is not a factor\n")
  globaln <- nrow(MicDtDF)
  if (length(agrby)!=globaln) stop("Size of the agrby argument does not agree with the number of rows in the MicDtDF data frame\n") 
#To  do: test if all columns of MicDtDF are numeric
#  if ( agrcrt[1]!="minmax" && (class(agrcrt)!="numeric" || length(agrcrt)!=2 || agrcrt[1]>=agrcrt[2] || agrcrt[1]<0. || agrcrt[2]>1.) )
  if ( agrcrt[1]!="minmax" && (class(agrcrt)[1]!="numeric" || length(agrcrt)!=2 || agrcrt[1]>=agrcrt[2] || agrcrt[1]<0. || agrcrt[2]>1.) )
    stop("Wrong value for the agrcrt argument\n( it should be either the string minmax or a two-dim vector\nof a prob. value for the lower percentile, followed by the prob. value for the upper percentile - \nex:c(0.05,0.95) )\n") 
   
  if (length(unique(agrby))!=length(levels(agrby)))  agrby <- factor(agrby)
  grplvls <- levels(agrby)
  lbDF <- ubDF <- data.frame(MicDtDF[1,])
  bndsDF <- cbind(lbDF,ubDF)

  ngrps <- length(grplvls)
  nvar <- ncol(MicDtDF)
  for (r in 1:ngrps) 
  { 
    grp <- grplvls[r]
    rind <- which(agrby==grp)
    for (c in 1:nvar) {
      if (agrcrt[1]=="minmax") {
        bndsDF[r,c] <- min(MicDtDF[rind,c])
        bndsDF[r,nvar+c] <- max(MicDtDF[rind,c])
      } else {
        bndsDF[r,c] <- quantile(MicDtDF[rind,c],probs=agrcrt[1])
        bndsDF[r,nvar+c] <- quantile(MicDtDF[rind,c],probs=agrcrt[2])
      }
    }  
  }
  res <- IData(bndsDF,Seq="AllLb_AllUb",VarNames=names(MicDtDF),ObsNames=grplvls) 
  DegInT <- which(apply(res@LogR,1,function(v) any(!is.finite(v))))
  nDegInT <- length(DegInT)
  if (nDegInT>0) {
    if (nDegInT==res@NObs) {
      warning("No Idata object was created because all units had some degenerate intervals")
      return(NULL)
    }
    if (nDegInT<10) {
      if (nDegInT==1) {
        wmsg <- paste("Data unit",res@ObsNames[DegInT],"was eliminated because it lead to some degenerate intervals")
      } else {
        wmsg <- paste(
          "Data units",paste(res@ObsNames[DegInT],collapse=", "),"were eliminated because they lead to some degenerate intervals",sep="\n"
        )
      }  
    } else {
      wmsg <- paste(nDegInT,"data units were eliminated because they lead to some degenerate intervals")
    }
    warning(wmsg)
    return(res[-DegInT,])
  }
  res
}

