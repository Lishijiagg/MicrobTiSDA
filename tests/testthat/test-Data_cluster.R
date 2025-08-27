test_that("Data.cluster", {
  data("OSLO.infant.data")
  data("OSLO.infant.meta")
  oslo_data_filt = Data.filter(Data = OSLO.infant.data,metadata = OSLO.infant.meta,
                               OTU_counts_filter_value = 5000, OTU_filter_value = 0.2,
                               Group_var = 'ID')

  oslo_data_intp = Data.interpolate(Data = oslo_data_filt,metadata = OSLO.infant.meta,
                                    Group_var = 'ID',
                                    Sample_ID = 'ID',
                                    Sample_Time = 'day')

  oslo_data_int = oslo_data_intp$Interpolated_Data
  oslo_meta_int = oslo_data_intp$Interpolated_Data_metadata

  oslo_data_trans = Data.trans(Data = oslo_data_int,
                               metadata = oslo_meta_int,
                               Group_var = 'Group')

  oslo_data_design = Design(metadata = oslo_meta_int,
                            Group_var = 'Group',
                            Sample_ID = 'ID',
                            Sample_Time = 'Time',
                            Pre_processed_Data = oslo_data_trans)

  oslo_model = Reg.SPLR(Data_for_Reg = oslo_data_design,
                        pre_processed_data = oslo_data_trans,
                        Knots = c(22,33),
                        unique_values = 10)

  oslo_model_pred = Pred.data(Fitted_models = oslo_model,
                              metadata = oslo_meta_int,
                              Group = 'Group',
                              Sample_Time = 'Time',
                              time_step = 1)

  oslo_model_clust = Data.cluster(predicted_data = oslo_model_pred,
                                  clust_method = 'average',
                                  dend_title_size = 12,
                                  font_size = 2.5)
  expect_s3_class(oslo_model_clust,"MicrobTiSDA.cluster")
  expect_true(length(oslo_model_clust) == 3)
})
