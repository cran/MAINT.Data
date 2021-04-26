context("AgrMcDt")

test_that("AgrMcDt (with min-max default) creates an IDAta object with the correct dimensions and attributes", {
  MicroDt <- data.frame(X=1:9,Y=9:1)
  agrfct <- factor(c("A","B","A","C","B","C","B","A","A"))
  Idt <- AgrMcDt(MicroDt,agrfct)

  expect_is(Idt,"IData")
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(NbMicroUnits(Idt),c(A=4,B=3,C=2))
  expect_equal(names(Idt),c("X","Y"))
} )

test_that("AgrMcDt performs a correct agregation into (min-max based) IData objects", {
  MicroDt <- data.frame(X=1:9,Y=9:1)
  agrfct <- factor(c("A","B","A","C","B","C","B","A","A"))
  Idt <- AgrMcDt(MicroDt,agrfct)
  TrueIDt <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),
                   VarNames=c("X","Y"),ObsNames = c("A","B","C"),
                   NbMicroUnits=c(A=4,B=3,C=2))
  
  expect_identical(Idt,TrueIDt)
} )

test_that("quantile based AgrMcDt creates an IDAta object with the correct dimensions and attributes", {
  MicroDt <- data.frame(X=rep(0:10,3),Y=rep(seq(0,100,by=10),3))
  agrfct <- factor(c(rep("A",6),rep("B",11),rep("A",5),rep("C",11)))
  Idt <- AgrMcDt(MicroDt,agrby=agrfct,agrcrt=c(0.1,0.9))

  expect_is(Idt,"IData")
  expect_equal(nrow(Idt),3)
  expect_equal(ncol(Idt),2)
  expect_equal(NbMicroUnits(Idt),c(A=11,B=11,C=11))
  expect_equal(names(Idt),c("X","Y"))
} )

test_that("AgrMcDt performs a correct agregation into quantile based IData objects", {
  MicroDt <- data.frame(X=rep(0:10,3),Y=rep(seq(0,100,by=10),3))
  agrfct <- factor(c(rep("A",6),rep("B",11),rep("A",5),rep("C",11)))
  Idt <- AgrMcDt(MicroDt,agrby=agrfct,agrcrt=c(0.1,0.9))
  TrueIDt <- IData(data.frame(list(rep(1,3),rep(9,3),rep(10,3),rep(90,3))),
                   VarNames=c("X","Y"),ObsNames = c("A","B","C"),
                   NbMicroUnits=c(A=11,B=11,C=11))
  
  expect_identical(Idt,TrueIDt)
} )

