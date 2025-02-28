## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----setup--------------------------------------------------------------------
#  library("MicrobTiSDA")
#  library("tidyr")
#  library("dplyr")
#  library("ggplot2")

## -----------------------------------------------------------------------------
#  data("OSLO.infant.data")
#  data("OSLO.infant.meta")
#  data("OSLO.infant.taxa")

## -----------------------------------------------------------------------------
#  oslo_data_filt = Data.filter(Data = OSLO.infant.data,metadata = OSLO.infant.meta,
#                               OTU_counts_filter_value = 5000, OTU_filter_value = 0.2,
#                               Group_var = 'ID')

## -----------------------------------------------------------------------------
#  oslo_data_intp = Data.interpolate(Data = oslo_data_filt,metadata = OSLO.infant.meta,
#                                    Group_var = 'ID',
#                                    Sample_ID = 'ID',
#                                    Sample_Time = 'day')

## -----------------------------------------------------------------------------
#  oslo_data_int = oslo_data_intp$Interpolated_Data
#  oslo_meta_int = oslo_data_intp$Interpolated_Data_metadata

## -----------------------------------------------------------------------------
#  oslo_data_trans = Data.trans(Data = oslo_data_int,
#                               metadata = oslo_meta_int,
#                               Group_var = 'Group')

## -----------------------------------------------------------------------------
#  oslo_data_interact = Spec.interact(Data = oslo_data_trans,
#                                     metadata = oslo_meta_int,
#                                     Group_var = 'Group',
#                                     num_iterations = 5)

## -----------------------------------------------------------------------------
#  oslo_interact_vis = Interact.vis(Interact_data = oslo_data_interact,
#                                   count_data = oslo_data_trans,
#                                   metadata = oslo_meta_int,
#                                   Subject_ID = 'ID',
#                                   Interact_threshold = 1e-6,
#                                   legend_text_size = 8,
#                                   legend_title_size = 10,
#                                   title_size = 12,
#                                   label_distance = 0.35,
#                                   label_font_size = 3)

## -----------------------------------------------------------------------------
#  oslo_data_design = Design(metadata = oslo_meta_int,
#                            Group_var = 'Group',
#                            Sample_ID = 'ID',
#                            Sample_Time = 'Time',
#                            Pre_processed_Data = oslo_data_trans)

## -----------------------------------------------------------------------------
#  oslo_model = Reg.SPLR(Data_for_Reg = oslo_data_design,
#                        pre_processed_data = oslo_data_trans,
#                        Knots = c(22,33),
#                        unique_values = 10)

## -----------------------------------------------------------------------------
#  oslo_model_pred = Pred.data(Fitted_models = oslo_model,
#                              metadata = oslo_meta_int,
#                              Group = 'Group',
#                              Sample_Time = 'Time',
#                              time_step = 1)

## -----------------------------------------------------------------------------
#  oslo_model_clust = Data.cluster(predicted_data = oslo_model_pred,
#                                  clust_method = 'average',
#                                  auto_cutree = FALSE,
#                                  dend_title_size = 12,
#                                  font_size = 2.5)

## ----message=TRUE-------------------------------------------------------------
#  oslo_model_vis = Data.visual(cluster_results = oslo_model_clust,
#                               cutree_by = 'height',
#                               cluster_height = c(0.15,0.15),
#                               predicted_data = oslo_model_pred,
#                               Design_data = oslo_data_design,
#                               pre_processed_data = oslo_data_trans,
#                               plot_dots = TRUE,
#                               figure_x_scale = 5,
#                               Taxa = OSLO.infant.taxa,
#                               plot_lm = FALSE,
#                               legend_title_size = 10,
#                               legend_text_size = 8)
#  
#  # Visualize the first OTU cluster in ID10
#  oslo_model_vis$ID_10[[1]]

## -----------------------------------------------------------------------------
#  opposite_vis = Data.opp.cor.vis(predicted_data = oslo_model_pred,
#                                  pre_processed_data = oslo_data_trans,
#                                  Design_data = oslo_data_design,
#                                  ng_cor_thres = 0.1,
#                                  Taxa = OSLO.infant.taxa,
#                                  plot_dots = FALSE,legend_text_size = 5)
#  
#  # For example, taking OTU_000001 of ID10
#  # we can observe microbial features with a correlation distance greater than 0.9
#  opposite_vis$ID_10$OTU_000001

## -----------------------------------------------------------------------------
#  data("preterm.data")
#  data("preterm.meta")
#  data("preterm.taxa")

## -----------------------------------------------------------------------------
#  preterm_rf = Data.rf.classifier(raw_data = preterm.data,
#                                  metadata = preterm.meta,
#                                  train_p = 0.7,
#                                  Group = 'Group',
#                                  OTU_counts_filter_value = 0)

## -----------------------------------------------------------------------------
#  preterm_rf_vis = Classify.vis(classified_results = preterm_rf,dist_method = 'bray')

## -----------------------------------------------------------------------------
#  preterm_trans = Data.trans(Data = preterm_rf$Important_OTU_table,
#                             Group_var = 'Group',
#                             metadata = preterm.meta)

## -----------------------------------------------------------------------------
#  preterm_design = Design(metadata = preterm.meta,
#                          Group_var = 'Group',
#                          Pre_processed_Data = preterm_trans,
#                          Sample_ID = 'ID',
#                          Sample_Time = 'Time')

## -----------------------------------------------------------------------------
#  preterm_model= Reg.MESR(Data_for_Reg = preterm_design,
#                          pre_processed_data = preterm_trans,
#                          max_Knots = 2)

## -----------------------------------------------------------------------------
#  preterm_pred = Pred.data.MESR(Fitted_models = preterm_model,
#                                metadata = preterm.meta,
#                                Group = 'Group',
#                                Sample_Time = 'Time',
#                                time_step = 1)

## -----------------------------------------------------------------------------
#  preterm_clust = Data.cluster(predicted_data = preterm_pred,
#                               clust_method = 'average',
#                               auto_cutree = FALSE,
#                               font_size = 3.5,
#                               dend_title_size = 12)

## -----------------------------------------------------------------------------
#  preterm_data_vis = Data.visual.MESR(Design_data = preterm_design,
#                                      cluster_results = preterm_clust,
#                                      cutree_by = 'height',
#                                      cluster_height = c(0.2,0.2),
#                                      predicted_data = preterm_pred,
#                                      pre_processed_data = preterm_trans,
#                                      Taxa = preterm.taxa,
#                                      plot_dots = TRUE,
#                                      plot_lm = FALSE,
#                                      figure_x_scale = 1,
#                                      legend_title_size = 10,
#                                      legend_text_size = 8,
#                                      axis_title_size = 10,
#                                      axis_x_size = 7,
#                                      axis_y_size = 7)
#  
#  # Example: visualize the first OTU cluster in both the S_C and S_D groups
#  preterm_data_vis$S_C[[1]]
#  preterm_data_vis$S_D[[1]]

