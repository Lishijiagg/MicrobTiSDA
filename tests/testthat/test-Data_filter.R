test_that("Data.filter", {
  set.seed(123)
  otu_table <- as.data.frame(matrix(sample(0:100, 100, replace = TRUE), nrow = 10))
  rownames(otu_table) <- paste0("OTU", 1:10)
  colnames(otu_table) <- paste0("Sample", 1:10)

  metadata <- data.frame(Group = rep(c("A", "B"), each = 5),row.names = paste0("Sample", 1:10))

  # 1. Test OTU_counts_filter_value
  result <- Data.filter(otu_table, metadata, OTU_counts_filter_value = 50, OTU_filter_value = 0.2)
  expect_true(all(rowSums(result$filtered_table) > 15))

  # 2. Test error report correctly
  expect_error(Data.filter(matrix(1:4, nrow = 2), test_metadata))

  # 3. Test Group_var
  result_group <- Data.filter(otu_table, metadata, OTU_counts_filter_value = 0, OTU_filter_value = 0.5, Group_var = "Group")
  expect_s3_class(result_group, "FilteredData")
})

