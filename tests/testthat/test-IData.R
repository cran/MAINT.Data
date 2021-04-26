context("IData")

test_that("Idata reads correctly a data frame of lower and upper bounds presented variable by variable", {

  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  Allbnds <- cbind(Xbnds,Ybnds)

  Idt <- IData(Allbnds,VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Allbnds))
  expect_equal(rownames(Idt),rownames(Allbnds))
  
  X.MidP <- (Xbnds$XLB + Xbnds$XUB) /2 
  Y.MidP <- (Ybnds$YLB + Ybnds$YUB) /2 
  AllMidP <- as.data.frame(cbind(X.MidP,Y.MidP))
  row.names(AllMidP) <- row.names(Allbnds)

  expect_identical(MidPoints(Idt),AllMidP)
    
  X.LogR <- log(Xbnds$XUB - Xbnds$XLB) 
  Y.LogR <- log(Ybnds$YUB - Ybnds$YLB) 
  AllLogR <- as.data.frame(cbind(X.LogR,Y.LogR))
  row.names(AllLogR) <- row.names(Allbnds)
  
  expect_identical(LogRanges(Idt),AllLogR)
  
} )

test_that("Idata reads correctly a data frame of all lower bounds followed by all upper bounds", {
  
  Lbnds <- data.frame(XLB=c(1,2,4),YLB=c(10,3,5),row.names=paste("Unit",1:3))
  Ubnds <- data.frame(XUB=c(9,7,6),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  Allbnds <- cbind(Lbnds,Ubnds)
  
  Idt <- IData(Allbnds,Seq="AllLb_AllUb",VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Allbnds))
  expect_equal(rownames(Idt),rownames(Allbnds))
  
  X.MidP <- (Lbnds$XLB + Ubnds$XUB) /2 
  Y.MidP <- (Lbnds$YLB + Ubnds$YUB) /2 
  AllMidP <- as.data.frame(cbind(X.MidP,Y.MidP))
  row.names(AllMidP) <- row.names(Allbnds)
  
  expect_identical(MidPoints(Idt),AllMidP)
  
  X.LogR <- log(Ubnds$XUB - Lbnds$XLB) 
  Y.LogR <- log(Ubnds$YUB - Lbnds$YLB) 
  AllLogR <- as.data.frame(cbind(X.LogR,Y.LogR))
  row.names(AllLogR) <- row.names(Allbnds)
  
  expect_identical(LogRanges(Idt),AllLogR)
  
} )

test_that("Idata reads correctly a data frame of MidPoints and LogRanges presented variable by variable", {
  
  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  X.MidP <- (Xbnds$XLB + Xbnds$XUB) /2 
  Y.MidP <- (Ybnds$YLB + Ybnds$YUB) /2 
  X.LogR <- log(Xbnds$XUB - Xbnds$XLB) 
  Y.LogR <- log(Ybnds$YUB - Ybnds$YLB) 
  MidLogRDF <- data.frame(cbind(X.MidP,X.LogR,Y.MidP,Y.LogR),row.names=rownames(Xbnds))
  
  Idt <- IData(MidLogRDF,Seq="MidPLogR_VarbyVar",VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Xbnds))
  expect_equal(rownames(Idt),rownames(Xbnds))
  
  AllMidP <- as.data.frame(cbind(X.MidP,Y.MidP))
  row.names(AllMidP) <- row.names(Xbnds)
  expect_identical(MidPoints(Idt),AllMidP)
  
  AllLogR <- as.data.frame(cbind(X.LogR,Y.LogR))
  row.names(AllLogR) <- row.names(Xbnds)
  expect_identical(LogRanges(Idt),AllLogR)
  
} )


test_that("Idata reads correctly a data frame of all MidPoints followed by all LogRanges", {
  
  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  X.MidP <- (Xbnds$XLB + Xbnds$XUB) /2 
  Y.MidP <- (Ybnds$YLB + Ybnds$YUB) /2 
  X.LogR <- log(Xbnds$XUB - Xbnds$XLB) 
  Y.LogR <- log(Ybnds$YUB - Ybnds$YLB) 

  AllMidP <- as.data.frame(cbind(X.MidP,Y.MidP))
  row.names(AllMidP) <- row.names(Xbnds)
  AllLogR <- as.data.frame(cbind(X.LogR,Y.LogR))
  row.names(AllLogR) <- row.names(Xbnds)
  MidLogRDF <- cbind(AllMidP,AllLogR)
  
  Idt <- IData(MidLogRDF,Seq="AllMidP_AllLogR",VarNames = c("X","Y"))
  
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(row.names(Idt),row.names(Xbnds))
  expect_equal(rownames(Idt),rownames(Xbnds))
  
  expect_identical(MidPoints(Idt),AllMidP)
  expect_identical(LogRanges(Idt),AllLogR)
  
} )

