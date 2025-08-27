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

  selected_biomarkers = Rf.biomarkers(rf = preterm_rf,feature_select_num = 10)
  expect_s3_class(selected_biomarkers,"RfBiomarker")
})
