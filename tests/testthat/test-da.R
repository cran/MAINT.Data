context("discriminant analysis methods")

test_that("lda works correctly", {

  Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
  AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                             ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
                             )
                  )

  ldares <- lda(AbaloneIdt,AbalClass)
  ldapred <- predict(ldares,AbaloneIdt)
  PPsums <- rowSums(ldapred$posterior)
  names(PPsums) <- NULL
  expect_equal( PPsums, rep(1.,nrow(AbaloneIdt)) )
  
  PPpred <- as.factor(levels(AbalClass)[apply(ldapred$posterior,1,which.max)])
  names(PPpred) <- rownames(AbaloneIdt)
  expect_equal(PPpred,ldapred$class)
  
  Trueldapred <- factor(c("Adult","Adult","Adult","Old","Old","Old","Young","Young","Young","Adult","Adult","Adult",
                          "Old","Young","Young","Young","Adult","Adult","Adult","Old","Old","Old","Young","Young"))
  names(Trueldapred) <- rownames(AbaloneIdt)
  expect_equal(Trueldapred,ldapred$class)
  
} )

test_that("qda works correctly", {
  
  Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
  AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                             ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
  )
  )
  
  qdares <- qda(AbaloneIdt[,2],AbalClass)        # Note: using lower-dimensionality data set in order to avoid singular covariance matrices
  qdapred <- predict(qdares,AbaloneIdt[,2])
  PPsums <- rowSums(qdapred$posterior)
  names(PPsums) <- NULL
  expect_equal( PPsums, rep(1.,nrow(AbaloneIdt)) )
  
  PPpred <- as.factor(levels(AbalClass)[apply(qdapred$posterior,1,which.max)])
  names(PPpred) <- rownames(AbaloneIdt[,2])
  expect_equal(PPpred,qdapred$class)
  
  Trueqdapred <- factor(c("Adult","Old","Old","Old","Old","Old","Young","Young","Young","Young","Young","Adult",
                          "Adult","Young","Young","Young","Adult","Adult","Old","Old","Old","Old","Young","Young"))
  names(Trueqdapred) <- rownames(AbaloneIdt)
  expect_equal(Trueqdapred,qdapred$class)
  
} )

test_that("leave-one-out cross-validation works correctly", {
   
   Agestrg <- substring(AbaloneIdt@ObsNames,first=3)
   AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"|Agestrg=="7-9","Young",
                              ifelse(Agestrg=="10-12"|Agestrg=="13-15"|Agestrg=="16-18","Adult","Old")
                             )
                      )
   
   looldares <- DACrossVal(AbaloneIdt,AbalClass,lda,loo=TRUE)
   Nk <- colSums(looldares[,,"Nk"])  
   TrueNk <- as.numeric(table(AbalClass))
   names(TrueNk) <- levels(AbalClass)
   expect_equal(Nk,TrueNk)

   looerrors <- colSums(looldares[,,"Nk"] * looldares[,,"Clerr"], na.rm=TRUE)  
   Turelooerrors <- c(2,3,4)
   names(Turelooerrors) <- levels(AbalClass)
   expect_equal(Turelooerrors,looerrors)
   
   looerrest <- looerrors /Nk
   Glooerrest <- sum(looerrors)/sum(Nk)
   names(Glooerrest) <- "Global"
   expect_equal(c(looerrest,Glooerrest),attr(looldares,"errestimates"))
   
} )
