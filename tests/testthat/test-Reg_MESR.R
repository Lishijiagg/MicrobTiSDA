test_that("Reg.MESR", {
  data <- data.frame(Feature1 = rnorm(120),
                     Feature2 = rnorm(120),
                     Feature3 = rnorm(120),
                     Feature4 = rnorm(120),
                     Feature5 = rnorm(120))

  rownames(data) = paste0("Sample",seq(1:120))

  metadata = data.frame(TimePoint = rep(1:20,times = 6),
                        Sample = paste0("Sample",seq(1:120)),
                        Group = rep(c("A","B"),each = 60))

  design_data = Design(metadata = metadata,
                       Group_var = 'Group',
                       Pre_processed_Data = data,
                       Sample_ID = 'Sample',
                       Sample_Time = 'TimePoint')
  design_data$data$ID <- as.factor(design_data$data$ID)

  fit_model = Reg.MESR(Data_for_Reg = design_data,
                       pre_processed_data = data,
                       max_Knots = 3,unique_values = 3)
  expect_s3_class(fit_model,"RegMESR")
  expect_true(length(fit_model$fitted_model) == length(unique(metadata$Group)))
  expect_true(length(fit_model$fitted_model$A) == ncol(data))
  expect_true(length(fit_model$fitted_model$B) == ncol(data))
})
