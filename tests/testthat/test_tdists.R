context("Fitting mixtures of t-distributions")

test_that("it at least runs", {
  data(alon_data)
  res<-each_gene(alon_data[1,])
  expect_equal(res$k, 3)
})

test_that("Does it at least run", {
    data(alon_data)
    set.seed(56)
    res<-each_gene(alon_data[2,])
    expect_equal(res$k, 3)
})