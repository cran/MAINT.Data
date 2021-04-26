context("binding functions")

test_that("cbind works correctly with two arguments", {
  Idt1 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"))
  Idt2 <- IData(data.frame(list(c(3,8,4),c(5,9,7))),VarNames=c("Z"))
  Idt12 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6),c(3,8,4),c(5,9,7))),VarNames=c("X","Y","Z"))
                                              
  expect_identical(cbind(Idt1,Idt2),Idt12)
} )

test_that("cbind works correctly with more than two arguments", {
  Idt1 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"))
  Idt2 <- IData(data.frame(list(c(3,8,4),c(5,9,7))),VarNames="Z")
  Idt3 <- IData(data.frame(list(c(2,5,6),c(4,10,8))),VarNames="W")
  Idt123 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6),c(3,8,4),c(5,9,7),c(2,5,6),c(4,10,8))),VarNames=c("X","Y","Z","W"))
  
  expect_identical(cbind(Idt1,Idt2,Idt3),Idt123)
} )


test_that("rbind works correctly with two arguments", {
  
  Idt1 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"), ObsNames = c("A","B","C"))
  Idt2 <- IData(data.frame(list(c(3,8,4),c(5,9,7)),c(2,6,3),c(7,8,4)),VarNames=c("X","Y"), ObsNames = c("D","E","F"))
  Idt12 <- IData(data.frame(list(c(1,2,4,3,8,4),c(9,7,6,5,9,7),c(1,3,4,2,6,3),c(9,8,6,7,8,4))),
                 VarNames=c("X","Y"),ObsNames = c("A","B","C","D","E","F"))
  
  expect_identical(rbind(Idt1,Idt2),Idt12)
  
  
  Idt1 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"))
  Idt2 <- IData(data.frame(list(c(3,8,4),c(5,9,7)),c(2,6,3),c(7,8,4)),VarNames=c("X","Y"))
  Idt12 <- IData(data.frame(list(c(1,2,4,3,8,4),c(9,7,6,5,9,7),c(1,3,4,2,6,3),c(9,8,6,7,8,4))),
                 VarNames=c("X","Y"),ObsNames = as.character(c(1,2,3,11,21,31)))
  
  expect_identical(rbind(Idt1,Idt2),Idt12)
} )

test_that("rbind works correctly with more than two arguments", {
  
  Idt1 <- IData(data.frame(list(c(1,2,4),c(9,7,6),c(1,3,4),c(9,8,6))),VarNames=c("X","Y"), ObsNames = c("A","B","C"))
  Idt2 <- IData(data.frame(list(c(3,8,4),c(5,9,7)),c(2,6,3),c(7,8,4)),VarNames=c("X","Y"), ObsNames = c("D","E","F"))
  Idt3 <- IData(data.frame(list(2,3,6,8)),VarNames=c("X","Y"), ObsNames = "G")
  Idt123 <- IData(data.frame(list(c(1,2,4,3,8,4,2),c(9,7,6,5,9,7,3),c(1,3,4,2,6,3,6),c(9,8,6,7,8,4,8))),
                 VarNames=c("X","Y"),ObsNames = c("A","B","C","D","E","F","G"))
  
  expect_identical(rbind(Idt1,Idt2,Idt3),Idt123)
  
} )
