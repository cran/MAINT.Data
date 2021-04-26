context("IData indexing operations")

test_that("row-only indexing works correctly for Idata objects", {

  Xbnds <- data.frame(XLB=c(1,2,4),XUB=c(9,7,6),row.names=paste("Unit",1:3))
  Ybnds <- data.frame(YLB=c(10,3,5),YUB=c(15,8,6),row.names=paste("Unit",1:3))
  Allbnds <- cbind(Xbnds,Ybnds)
  Idt <- IData(Allbnds,VarNames = c("X","Y"))
  
  Idt2 <- IData(Allbnds[2,],VarNames = c("X","Y"))
  Idt13 <- IData(Allbnds[c(1,3),],VarNames = c("X","Y"))
  Idt31 <- IData(Allbnds[c(3,1),],VarNames = c("X","Y"))
  Idt23 <- IData(Allbnds[2:3,],VarNames = c("X","Y"))
  Idt32 <- IData(Allbnds[3:2,],VarNames = c("X","Y"))
  
  expect_identical(Idt[2,],Idt2)
  expect_identical(Idt[-c(1,3),],Idt2)
  expect_identical(Idt[c(1,3),],Idt13)
  expect_identical(Idt[-2,],Idt13)
  expect_identical(Idt[c(3,1),],Idt31)
  expect_identical(Idt[2:3,],Idt23)
  expect_identical(Idt[-1,],Idt23)
  expect_identical(Idt[3:2,],Idt32)

} )


test_that("column-only indexing works correctly for Idata objects", {
  
  Lbnds <- data.frame(XLB=c(1,2,4),YLB=c(10,3,5),ZLB=c(8,5,4),row.names=paste("Unit",1:3))
  Ubnds <- data.frame(XUB=c(9,7,6),YUB=c(15,8,6),ZUB=c(12,9,7),row.names=paste("Unit",1:3))
  Idt <- IData(cbind(Lbnds,Ubnds),Seq="AllLb_AllUb",VarNames = c("X","Y","Z"))
  
  Idt2 <- IData(cbind(Lbnds[,2,drop=FALSE],Ubnds[,2,drop=FALSE]),Seq="AllLb_AllUb",VarNames = "Y")
  Idt13 <- IData(cbind(Lbnds[,c(1,3)],Ubnds[,c(1,3)]),Seq="AllLb_AllUb",VarNames = c("X","Z"))
  Idt31 <- IData(cbind(Lbnds[,c(3,1)],Ubnds[,c(3,1)]),Seq="AllLb_AllUb",VarNames = c("Z","X"))
  Idt23 <- IData(cbind(Lbnds[,2:3],Ubnds[,2:3]),Seq="AllLb_AllUb",VarNames = c("Y","Z"))
  Idt32 <- IData(cbind(Lbnds[,3:2],Ubnds[,3:2]),Seq="AllLb_AllUb",VarNames = c("Z","Y"))
  
  expect_identical(Idt[,2],Idt2)
  expect_identical(Idt[,-c(1,3)],Idt2)
  expect_identical(Idt[,c(1,3)],Idt13)
  expect_identical(Idt[,-2],Idt13)
  expect_identical(Idt[,c(3,1)],Idt31)
  expect_identical(Idt[,2:3],Idt23)
  expect_identical(Idt[,-1],Idt23)
  expect_identical(Idt[,3:2],Idt32)
  
} )


test_that("indexing by both rows and columns works correctly for Idata objects", {
  
  Lbnds <- data.frame(XLB=c(1,2,4),YLB=c(10,3,5),ZLB=c(8,5,4),row.names=paste("Unit",1:3))
  Ubnds <- data.frame(XUB=c(9,7,6),YUB=c(15,8,6),ZUB=c(12,9,7),row.names=paste("Unit",1:3))
  Idt <- IData(cbind(Lbnds,Ubnds),Seq="AllLb_AllUb",VarNames = c("X","Y","Z"))
  
  Idtr2c3 <- IData(cbind(Lbnds[2,3,drop=FALSE],Ubnds[2,3,drop=FALSE]),Seq="AllLb_AllUb",VarNames = "Z")
  Idtr2c13 <- IData(cbind(Lbnds[2,c(1,3),drop=FALSE],Ubnds[2,c(1,3)]),Seq="AllLb_AllUb",VarNames = c("X","Z"))
  Idtr13c2 <- IData(cbind(Lbnds[c(1,3),2,drop=FALSE],Ubnds[c(1,3),2]),Seq="AllLb_AllUb",VarNames = "Y")
  Idtr13c21 <- IData(cbind(Lbnds[c(1,3),2:1],Ubnds[c(1,3),2:1]),Seq="AllLb_AllUb",VarNames = c("Y","X"))
  Idtr32c12 <- IData(cbind(Lbnds[3:2,1:2],Ubnds[3:2,1:2]),Seq="AllLb_AllUb",VarNames = c("X","Y"))
  
  expect_identical(Idt[2,3],Idtr2c3)
  expect_identical(Idt[-c(1,3),3],Idtr2c3)
  expect_identical(Idt[2,-(1:2)],Idtr2c3)
  expect_identical(Idt[-c(1,3),-c(1,2)],Idtr2c3)

  expect_identical(Idt[2,c(1,3)],Idtr2c13)
  expect_identical(Idt[-c(1,3),c(1,3)],Idtr2c13)
  expect_identical(Idt[2,-2],Idtr2c13)
  expect_identical(Idt[-c(1,3),-2],Idtr2c13)
  
  expect_identical(Idt[c(1,3),2],Idtr13c2)
  expect_identical(Idt[-2,2],Idtr13c2)
  expect_identical(Idt[c(1,3),-c(1,3)],Idtr13c2)
  expect_identical(Idt[-2,-c(1,3)],Idtr13c2)
  
  expect_identical(Idt[c(1,3),2:1],Idtr13c21)
  expect_identical(Idt[-2,2:1],Idtr13c21)

  expect_identical(Idt[3:2,1:2],Idtr32c12)
  expect_identical(Idt[3:2,-3],Idtr32c12)

} )
