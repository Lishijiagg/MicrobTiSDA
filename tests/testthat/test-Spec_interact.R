test_that("Spec.interact", {
  set.seed(123)
  Data <- as.data.frame(matrix(sample(1:100, 100, replace = TRUE), nrow = 10))
  colnames(Data) <- paste0("Feature", 1:10)
  rownames(Data) <- paste0("Sample", 1:10)

  metadata <- data.frame(Group = rep(c("A", "B"), each = 5),
                         Time = rep(c(1,2,3,4,5,1,2,3,4,5)))
  rownames(metadata) <- paste0("Sample", 1:10)

  results <- Spec.interact(Data = Data,
                           metadata = metadata,
                           Group_var = "Group",
                           abund_centered_method = "median",
                           num_iterations = 2,
                           error_threshold = 1e-3,
                           pre_error = 10000)

  expect_true(length(results$results) == length(unique(metadata$Group)))
  expect_error(Spec.interact(Data = Data,metadata = metadata,Group_var = "Group",abund_centered_method = NA,
                             num_iterations = 2,error_threshold = 1e-3,pre_error = 10000))
  expect_s3_class(results,"microbTiSDA")

})
