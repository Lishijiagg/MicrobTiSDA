test_that("multiplication works", {
  data("preterm.data")
  data("preterm.meta")
  preterm_rf = Data.rf.classifier(raw_data = preterm.data,
                                  metadata = preterm.meta,
                                  train_p = 0.7,
                                  Group = 'Group',
                                  OTU_counts_filter_value = 0,
                                  legend_title_size = 12,
                                  legend_text_size = 10,
                                  axis_title_size = 10,
                                  title_size = 12)
  expect_s3_class(preterm_rf,"DataRFClassifier")
  expect_true(dim(preterm_rf$Input_data)[1] == dim(preterm.data)[2])
  expect_true((dim(preterm_rf$Input_data)[2]-1) == dim(preterm.data)[1])
  expect_true(nrow(preterm.data) == nrow(preterm_rf$OTU_importance))
})
