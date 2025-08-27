test_that("multiplication works", {
  data("OSLO.infant.data")
  data("OSLO.infant.meta")
  data("OSLO.infant.taxa")
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

  opposite_vis = Data.opp.cor.vis(predicted_data = oslo_model_pred,
                                  pre_processed_data = oslo_data_trans,
                                  Design_data = oslo_data_design,
                                  ng_cor_thres = 0.1,
                                  Taxa = OSLO.infant.taxa,
                                  plot_dots = FALSE,
                                  legend_text_size = 5)
  expect_s3_class(opposite_vis,"DataOppCorVis")
  expect_true(length(opposite_vis$curves_plot$ID_11) == length(oslo_model$fitted_model$ID_11))
  expect_true(length(opposite_vis$curves_plot$ID_10) == length(oslo_model$fitted_model$ID_10))
})
