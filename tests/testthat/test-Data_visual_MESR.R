test_that("Data.visual.MESR", {
  data <- data.frame(Feature1 = rnorm(120),
                     Feature2 = rnorm(120),
                     Feature3 = rnorm(120),
                     Feature4 = rnorm(120),
                     Feature5 = rnorm(120))

  rownames(data) <- paste0("Sample",seq(1:120))

  metadata <- data.frame(TimePoint = rep(1:20,times = 6),
                         ID = paste0("Sample",seq(1:120)),
                         Group = rep(c("A","B"),each = 60))

  design_data <- Design(metadata = metadata,
                        Group_var = 'Group',
                        Pre_processed_Data = data,
                        Sample_ID = 'ID',
                        Sample_Time = 'TimePoint')
  design_data$data$ID <- as.factor(design_data$data$ID)

  fit_model <- Reg.MESR(Data_for_Reg = design_data,
                        pre_processed_data = data,
                        max_Knots = 3,unique_values = 3)
  model_pred <- Pred.data.MESR(Fitted_models = fit_model,
                               metadata = metadata,
                               Group = 'Group',
                               Sample_Time = 'TimePoint',
                               time_step = 1)
  clust <- Data.cluster(predicted_data = model_pred,
                        clust_method = 'average',
                        dend_title_size = 12,
                        font_size = 3.5)

  clust_cut <- Data.cluster.cut(cluster_outputs = clust,cut_height = 0.15,font_size = 3.5)

  clust_plot <- Data.visual.MESR(Design_data = design_data,
                                 cluster_results = clust_cut,
                                 cutree_by = 'height',
                                 cluster_height = c(0.15,0.15),
                                 predicted_data = model_pred,
                                 pre_processed_data = data,
                                 Taxa = NULL,
                                 plot_dots = TRUE,
                                 plot_lm = FALSE,
                                 figure_x_scale = 1,
                                 legend_title_size = 10,
                                 legend_text_size = 8,
                                 axis_title_size = 10,
                                 axis_x_size = 7,
                                 axis_y_size = 7)

  expect_true(length(clust_plot$plots) == length(unique(metadata$Group)))

  expect_error(Data.visual.MESR(Design_data = design_data,
                                cluster_results = clust_cut,
                                cutree_by = NA,
                                cluster_height = c(0.15,0.15),
                                predicted_data = model_pred,
                                pre_processed_data = data,
                                Taxa = NULL,
                                plot_dots = TRUE,
                                plot_lm = FALSE,
                                figure_x_scale = 1,
                                legend_title_size = 10,
                                legend_text_size = 8,
                                axis_title_size = 10,
                                axis_x_size = 7,
                                axis_y_size = 7))
})
