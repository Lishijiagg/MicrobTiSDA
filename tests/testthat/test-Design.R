test_that("Design", {
  metadata <- data.frame(TimePoint = c(1, 2, 3, 4),
                         Sample = c('S1', 'S2', 'S3', 'S4'),
                         GroupA = c('A', 'A', 'B', 'B'),
                         GroupB = c('X', 'Y', 'X', 'Y'))

  Pre_processed_Data <- data.frame(Feature1 = rnorm(4),
                                   Feature2 = rnorm(4))

  design_data <- Design(metadata,
                        Group_var = c('GroupA', 'GroupB'),
                        Pre_processed_Data,
                        Sample_Time = 'TimePoint',
                        Sample_ID = 'Sample')
  expect_s3_class(design_data,"Design")

  design_data_no_group <- Design(metadata,
                                 Group_var = NULL,
                                 Pre_processed_Data,
                                 Sample_Time = 'TimePoint',
                                 Sample_ID = 'Sample')
  expect_s3_class(design_data_no_group,"Design")
})
