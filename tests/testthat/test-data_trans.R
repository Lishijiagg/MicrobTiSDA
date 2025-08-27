test_that("Data.trans", {
  set.seed(123)
  Data <- matrix(sample(1:100, 50, replace = TRUE), nrow = 5)
  rownames(Data) <- paste0("Feature", 1:5)
  colnames(Data) <- paste0("Sample", 1:10)

  metadata <- data.frame(Group = rep(c("A", "B"), each = 5))
  rownames(metadata) <- paste0("Sample", 1:10)

  # Apply MCLR transformation to the entire dataset
  transformed_data <- Data.trans(Data, metadata, Group_var = NULL)
  expect_s3_class(transformed_data,"TransformedData")

  # Apply MCLR transformation separately for each group
  transformed_data_by_group <- Data.trans(Data, metadata, Group_var = "Group")
  expect_s3_class(transformed_data_by_group,"TransformedData")
})
