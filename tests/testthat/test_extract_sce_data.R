context('.extract_sce_data')

library(ggsc)
sce <- scuttle::mockSCE()
genes <- rownames(sce) |> sample(6)
samples <- colnames(sce) |> sample(6)

da <- ggsc:::.extract_sce_data(sce, features=genes, dims=NULL)

test_that("the specified features will be extract",{
  flag <- genes %in% colnames(da)
  testthat::expect_true(all(flag))  
})


da2 <- ggsc:::.extract_sce_data(sce, features = genes, 
                                dims = NULL, cells = samples)

test_that('The specified cells will be extract',{
  flag1 <- length(samples) == nrow(da2)
  flag2 <- samples %in% rownames(da2)
  testthat::expect_true(all(flag1, flag2))
})
