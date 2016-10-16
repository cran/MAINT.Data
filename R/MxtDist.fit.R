resdlmfit <- function(Data,grouping,freq) 
{
  if (!is.matrix(Data)) Data <- as.matrix(Data) 
  levels(grouping) <- levels(unique(grouping))
  X <- model.matrix(Data ~ grouping)
  lm.fit(X, Data, method="qr")  
}

MxtDist.fit <- function(EstFunc,Data,grouping,freq=rep(1,nrow(Data)),...) 
{
  X0 <- model.matrix(rep(1,length(levels(grouping))) ~ factor(levels(grouping)))
  fit0 <- resdlmfit(Data,grouping,freq)
  result  <- EstFunc(matrix(resid(fit0),nrow(Data),ncol(Data)),...)
  result$mures <- result$mu
  result$ksires <- result$ksi
  result$mu0 <- X0 %*% as.matrix(coef(fit0))
  result$mu <- scale(result$mu0,center=-result$mu,scale=FALSE)
  result$ksi <- scale(result$mu0,center=-result$ksi,scale=FALSE)
  attr(result$mu,"scaled:center") <- attr(result$ksi,"scaled:center") <- NULL
  result
}



