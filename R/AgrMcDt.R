AgrMcDt <- function(MicDtDF,agrby,agrcrt="minmax")
{
  if (!(is.data.frame(MicDtDF))) stop("First argument of AgMicroData must be a data frame\n")
  if (class(agrby)!="factor") stop("Argument agrby is not a factor\n")
  globaln <- nrow(MicDtDF)
  if (length(agrby)!=globaln) stop("Size of the agrby argument does not agree with the number of rows in the MicDtDF data frame\n") 
#To  do: test if all columns of MicDtDF are numeric
  if ( agrcrt[1]!="minmax" && (class(agrcrt)!="numeric" || length(agrcrt)!=2 || agrcrt[1]>=agrcrt[2] || agrcrt[1]<0. || agrcrt[2]>1.) )
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
  IData(bndsDF,Seq="AllLb_AllUb",VarNames=names(MicDtDF),ObsNames=grplvls) 
}

