context("Single Gene plots")

test_that("Plots work", {
    data(alon_data)
    set.seed(123)
    expect_equal_to_reference(plot_single_gene(alon_data,474), "test_1.rds") 
    expect_equal_to_reference(plot_single_gene(alon_data,1758, g=2),"test_2.rds") 
    expect_equal_to_reference(plot_single_gene(alon_data,1758,g=3),"test_3.rds") 
})
