EvalClrule <- function(darule,VData,Vgrp,k)
{
  grpcodes <- levels(Vgrp)
  errates <- array(dim=k)
  Nk <- array(dim=k)
  clres <- predict(darule,VData)$class
  if (k>2)   Cmat <- matrix(nrow=k,ncol=k,dimnames=list(grpcodes,grpcodes))
  else Cmat <- NULL
  for (grpInd in 1:k)
  {
    #    Nk[grpInd] <- length(Vgrp[Vgrp==grpcodes[grpInd]])
    Nk[grpInd] <- length(which(Vgrp==grpcodes[grpInd]))
    thisgrpclres <- clres[Vgrp==grpcodes[grpInd]]
    levels(thisgrpclres) <- grpcodes
    thisgrperr <- thisgrpclres[grpcodes[grpInd]!=thisgrpclres]
    if (Nk[grpInd]>0) {
      errates[grpInd] <- length(thisgrperr)/Nk[grpInd]
    } else {
      errates[grpInd] <- 0
    }
    for (grpInd2 in 1:k) Cmat[grpInd,grpInd2] <- length(which(thisgrpclres==grpcodes[grpInd2]))
  }
  list(err=errates,Nk=Nk,Cmat=Cmat)  #  return(list(err=errates,Nk=Nk,Cmat=Cmat))
}
