context("Fitting mixtures of t-distributions")

test_that("Does it at least run", {
  data(alon_data)
  res<-each_gene(alon_data[1,])
  expect_equal(res$k, 3)
  
})