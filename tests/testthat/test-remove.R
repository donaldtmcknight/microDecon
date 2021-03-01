test_that("load dataset", {
  
  load( system.file("data/Example_1.RData", package = "microDecon"))
  
  expect_type(Example_1,"list")
  expect_equal(dim(Example_1),c(6,8))
  
  
})


#run checks using devtools::test()