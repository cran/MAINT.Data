context("MANOVA")

test_that("MANOVA works correctly for default homoscedastic Gaussian models", {

  Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
  AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                             ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
                             )
                  )

  n <- nrow(AbaloneIdt)
  q <- 2*ncol(AbaloneIdt)            # Total number of MidPoints and LogRanges
  k <- length(levels(AbalClass))
    
  for (Cv in 1:4) {
    AbMANOVA <- MANOVA(AbaloneIdt,AbalClass,CovCase=Cv)
    AbMANOVAH0 <- H0res(AbMANOVA)
    AbaIdtmle <- mle(AbaloneIdt,CovCase=Cv)
    expect_equal(AbaIdtmle,AbMANOVAH0)
    AbMANOVAH1 <- H1res(AbMANOVA)
    SSCP <- matrix(0.,nrow=q,ncol=q)
    rownames(SSCP) <- colnames(SSCP) <- c(names(MidPoints(AbaloneIdt)),names(LogRanges(AbaloneIdt)))
    for (grp in levels(AbalClass)) {
      grpIdt <- AbaloneIdt[AbalClass==grp,]
      grpmle <- mle(grpIdt,CovCase=Cv)
      expect_equal(mean(grpmle),mean(AbMANOVAH1)[grp,])
      SSCP <- SSCP + nrow(grpIdt)*var(grpmle)
    }
    S <- SSCP/n
    expect_equal(S,var(AbMANOVAH1))
  }
  
} )

test_that("MANOVA  computes correct standar errors for default homoscedastic Gaussian models", {
  
  Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
  AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                             ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
  )
  )
  
  n <- nrow(AbaloneIdt)
  q <- 2*ncol(AbaloneIdt)            # Total number of MidPoints and LogRanges
  k <- length(levels(AbalClass))
  
  for (Cv in 1:4) {
    AbMANOVA <- MANOVA(AbaloneIdt,AbalClass,CovCase=Cv)
    AbMANOVAH0 <- H0res(AbMANOVA)
    
    AbmeanStder <- sd(AbMANOVAH0) / sqrt(n)
    expect_equal(stdEr(AbMANOVAH0)$mu,AbmeanStder)
    vcovb_AbmeanStder <- sqrt(diag(vcov(AbMANOVAH0)[1:q,1:q]))
    names(vcovb_AbmeanStder) <- names(AbmeanStder) <- NULL
    expect_equal(vcovb_AbmeanStder,AbmeanStder)
  
    mlecov <- var(AbMANOVAH0)
    mlecov[mlecov==0.] <- NA
    mlevar <- diag(mlecov)
    AbcovStder <- sqrt( (mlecov^2 + outer(mlevar,mlevar)) / (n-1) )   
    expect_equal(stdEr(AbMANOVAH0)$Sigma,AbcovStder)
  
    if (Cv==1) {  # Implement later vcov checks for other configurations
      vcovb_AbcovStder <- matrix(nrow=q,ncol=q)
      vcovb_AbcovStder[lower.tri(vcovb_AbcovStder,diag=TRUE)] <- sqrt(diag(vcov(AbMANOVAH0)[-(1:q),-(1:q)]))
      vcovb_AbcovStder[upper.tri(vcovb_AbcovStder)] <- t(vcovb_AbcovStder)[upper.tri(t(vcovb_AbcovStder))]
      dimnames(vcovb_AbcovStder) <- dimnames(AbcovStder)
      expect_equal(vcovb_AbcovStder,AbcovStder)
    }  
    
    AbMANOVAH1 <- H1res(AbMANOVA)
    nk <- as.numeric(table(AbalClass))

    SSCP <- matrix(0.,nrow=q,ncol=q)
    rownames(SSCP) <- colnames(SSCP) <- c(names(MidPoints(AbaloneIdt)),names(LogRanges(AbaloneIdt)))
    for (grp in levels(AbalClass)) {
      grpIdt <- AbaloneIdt[AbalClass==grp,]
      grpmle <- mle(grpIdt,CovCase=Cv)
      SSCP <- SSCP + nrow(grpIdt)*var(grpmle)
    }
    S <- SSCP/n
    AbH1meannames <- list( levels(AbalClass), c(names(MidPoints(AbaloneIdt)),names(LogRanges(AbaloneIdt))) )
    AbmeanStder <- matrix( rep(sqrt(diag(S)),each=k) / rep(sqrt(nk),q), nrow=k, ncol=q, dimnames=AbH1meannames )
    expect_equal(AbmeanStder,stdEr(AbMANOVAH1)$mu)
    vcovb_AbmeanStder <- matrix( sqrt(diag(vcov(AbMANOVAH1)[1:(k*q),1:(k*q)])), byrow=TRUE, nrow=k, ncol=q, dimnames=AbH1meannames )
    expect_equal(vcovb_AbmeanStder,AbmeanStder)
     
    mlecov <- var(AbMANOVAH1)
    mlecov[mlecov==0.] <- NA
    mlevar <- diag(mlecov)
    AbcovStder <- sqrt( (mlecov^2 + outer(mlevar,mlevar)) / (n-k) )   
    expect_equal(stdEr(AbMANOVAH1)$Sigma,AbcovStder)
     
    if (Cv==1) {  # Implement later vcov checks for other configurations
      vcovb_AbcovStder <- matrix(nrow=q,ncol=q)
      vcovb_AbcovStder[lower.tri(vcovb_AbcovStder,diag=TRUE)] <- sqrt(diag(vcov(AbMANOVAH1)[-(1:(k*q)),-(1:(k*q))]))
      vcovb_AbcovStder[upper.tri(vcovb_AbcovStder)] <- t(vcovb_AbcovStder)[upper.tri(t(vcovb_AbcovStder))]
      dimnames(vcovb_AbcovStder) <- dimnames(AbcovStder)
      expect_equal(vcovb_AbcovStder,AbcovStder)
    }  
  }
  
} )  
  
test_that("MANOVA works correctly for heteroscedastic Gaussian models", {
    
    Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
    AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                               ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
                        )
    )
    
    n <- nrow(AbaloneIdt)
    q <- 2*ncol(AbaloneIdt)            # Total number of MidPoints and LogRanges
    k <- length(levels(AbalClass))
    
    for (Cv in 1:4) {
      AbMANOVA <- MANOVA(AbaloneIdt,AbalClass,Mxt="Het",CovCase=Cv)
      AbMANOVAH0 <- H0res(AbMANOVA)
      AbaIdtmle <- mle(AbaloneIdt,CovCase=Cv)
      expect_equal(AbaIdtmle,AbMANOVAH0)
      AbMANOVAH1 <- H1res(AbMANOVA)
      for (grp in levels(AbalClass)) {
        grpIdt <- AbaloneIdt[AbalClass==grp,]
        grpmle <- mle(grpIdt,CovCase=Cv)
        expect_equal(mean(grpmle),mean(AbMANOVAH1)[grp,])
        expect_equal(var(grpmle),var(AbMANOVAH1)[,,grp])
      }
    }  
    
} )
  
test_that("MANOVA  computes correct standar errors for heteroscedastic Gaussian models", {
    
  Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
  AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                              ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
                            )
                     )
    
    n <- nrow(AbaloneIdt)
    q <- 2*ncol(AbaloneIdt)            # Total number of MidPoints and LogRanges
    k <- length(levels(AbalClass))
    
    for (Cv in 1:4) {
      AbMANOVA <- MANOVA(AbaloneIdt,AbalClass,Mxt="Het",CovCase=Cv)
      AbMANOVAH0 <- H0res(AbMANOVA)
      
      AbmeanStder <- sd(AbMANOVAH0) / sqrt(n)
      expect_equal(stdEr(AbMANOVAH0)$mu,AbmeanStder)
      vcovb_AbmeanStder <- sqrt(diag(vcov(AbMANOVAH0)[1:q,1:q]))
      names(vcovb_AbmeanStder) <- names(AbmeanStder) <- NULL
      expect_equal(vcovb_AbmeanStder,AbmeanStder)
      
      mlecov <- var(AbMANOVAH0)
      mlecov[mlecov==0.] <- NA
      mlevar <- diag(mlecov)
      AbcovStder <- sqrt( (mlecov^2 + outer(mlevar,mlevar)) / (n-1) )   
      expect_equal(stdEr(AbMANOVAH0)$Sigma,AbcovStder)
      
      if (Cv==1) {  # Implement later vcov checks for other configurations
        vcovb_AbcovStder <- matrix(nrow=q,ncol=q)
        vcovb_AbcovStder[lower.tri(vcovb_AbcovStder,diag=TRUE)] <- sqrt(diag(vcov(AbMANOVAH0)[-(1:q),-(1:q)]))
        vcovb_AbcovStder[upper.tri(vcovb_AbcovStder)] <- t(vcovb_AbcovStder)[upper.tri(t(vcovb_AbcovStder))]
        dimnames(vcovb_AbcovStder) <- dimnames(AbcovStder)
        expect_equal(vcovb_AbcovStder,AbcovStder)
      }  
      
      AbMANOVAH1 <- H1res(AbMANOVA)
      nk <- as.numeric(table(AbalClass))
      
      AbH1meannames <- list( levels(AbalClass), c(names(MidPoints(AbaloneIdt)),names(LogRanges(AbaloneIdt))) )
      AbmeanStder <- matrix( nrow=k, ncol=q, dimnames=AbH1meannames )
      for (g in 1:k) AbmeanStder[g,] <- sqrt( diag(var(AbMANOVAH1)[,,g]) / nk[g] ) 
      expect_equal(AbmeanStder,stdEr(AbMANOVAH1)$mu)
      vcovb_AbmeanStder <- matrix( nrow=k, ncol=q, dimnames=AbH1meannames )
      for (g in 1:k) vcovb_AbmeanStder[g,] <- sqrt( diag(suppressWarnings(vcov(AbMANOVAH1))[1:q,1:q,g]) ) 
      expect_equal(vcovb_AbmeanStder,AbmeanStder)
      
      mlecov <- var(AbMANOVAH1)
      mlecov[mlecov==0.] <- NA
      for (g in 1:k) {
        mlevar <- diag(mlecov[,,g])
        AbcovStder <- sqrt( (mlecov[,,g]^2 + outer(mlevar,mlevar)) / (nk[g]-1) )   
        expect_equal(stdEr(AbMANOVAH1)$Sigma[,,g],AbcovStder)
        if (Cv==1) {  # Implement later vcov checks for other configurations
          vcovb_AbcovStder <- matrix(nrow=q,ncol=q)
          vcovb_AbcovStder[lower.tri(vcovb_AbcovStder,diag=TRUE)] <- sqrt( diag(suppressWarnings(vcov(AbMANOVAH1))[-(1:q),-(1:q),g]) )
          vcovb_AbcovStder[upper.tri(vcovb_AbcovStder)] <- t(vcovb_AbcovStder)[upper.tri(t(vcovb_AbcovStder))]
          dimnames(vcovb_AbcovStder) <- dimnames(AbcovStder)
          expect_equal(vcovb_AbcovStder,AbcovStder)
        }  
      }  
      
    }  
    
} )
