test_that("Data.interpolate", {
  set.seed(123)
  Data <- matrix(sample(1:100, 40, replace = TRUE), nrow = 5)
  rownames(Data) <- paste0("Feature", 1:5)
  colnames(Data) <- paste0("Sample", 1:8)

  metadata <- data.frame(Time = c(1, 3, 5, 7, 2, 4, 6, 8),
                         ID = c(rep("Subject1", 4),
                                rep("Subject2", 4)),
                         Group = c(rep("A", 4),
                                   rep("B", 4)),
                         row.names = paste0("Sample", 1:8))
  interp_results <- Data.interpolate(Data = Data,
                                     metadata = metadata,
                                     Sample_Time = "Time",
                                     Sample_ID = "ID",
                                     interp_method = "cubic",
                                     Group_var = "Group")

  expect_s3_class(interp_results,"MicrobTiSDA.interpolate")
})
