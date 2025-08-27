test_that("Pred.data", {
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

  result <- Reg.SPLR(Data_for_Reg = design_data,
                     pre_processed_data = Pre_processed_Data,
                     z_score = NA,
                     unique_values = 0,
                     Knots = NULL,
                     max_Knots = 1)

  predictions <- Pred.data(result,metadata,
                           Group = "GroupA",
                           time_step = 1,
                           Sample_Time = "TimePoint")

  expect_true(length(predictions) == 0)
})
