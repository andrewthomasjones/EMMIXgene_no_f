context("Selecting Genes")

test_that("Number selected is consistent", {
    set.seed(123)
    alon_sel<-select_genes(alon_data[1:100,])
    n_alon<-sum(alon_sel$selected)
    expect_equal(n_alon, 12)
})
    
