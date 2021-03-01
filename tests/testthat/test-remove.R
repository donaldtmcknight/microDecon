test_that("load dataset", {
  
  load( sprintf("../../data/Example_1.Rdata",getwd()))
  
  expect_type(Example_1,"list")
  expect_equal(dim(Example_1),c(6,8))
  
  
})
