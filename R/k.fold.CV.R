#' @title Compare the Performance of MicrobTiSDA with Several Other Regression Models
#' @description
#' This function performs k-fold cross-validation on different regression models and compare their performance in
#'     predicting microbial feature time series data using root mean square error (RMSE).
#'
#' @details
#' The function randomly splits the dataset into k parts based on the k-fold parameter specified by the user. For
#'     each fold, the dataset is divided into training and testing subsets, making sure that the time range of the
#'     test set falls within the time range of the training set. Four regression models are then applied for fitting
#'     and prediction: (1) Polynomial regression fitting the relationship between time and microbial feature using the
#'     \code{MasigPro} package, with stepwise selection to determine the best model. (2) Using the \code{MetaDprof}
#'     package's ssanova method to fit a smoothing spline model. (3) A locally weighted regression model fitted using
#'     the \code{MetaLonDA} package. (4) A model fitted by trying degrees of freedom from 1 to 5 and choosing the optimal
#'     degree based on GCV, followed by prediction. For each model, the function calculates predictions on the test set
#'     and computes the RMSE by comparing the predictions with the actural values. he RMSE values are then compared to
#'     evaluate the accuracy of the different methods in predicting time-series microbial feature data.
#'
#'
#' @param Data_for_Reg The data frame output from the \code{\link[MicrobTiSDA]{Design}}.
#' @param pre_processed_data The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
#'     pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.
#' @param unique_values An integer specifying the minimum number of unique values required in an OTU for spline fitting (default: 5).
#' @param poly_degree A numeric value indicating the degree of the polynomial used in the GLM model (default: 3).
#' @param k An integer specifying the number of folds for cross-validation (default: 10).
#' @param box_plot Logical; if \code{TRUE}, boxplots of the CV RMSE values are generated (default: \code{TRUE}).
#' @param sig_label Logical; if \code{TRUE} (and if \code{box_plot = TRUE}), statistical significance annotations are
#'     added to the boxplots (default: \code{FALSE}).
#'
#' @return If \code{box_plot = TRUE}, a list with two elements: an RMSE table summarizing the mean RMSE of each method for each OTU,
#'     and a list of ggplot2 boxplots of the CV RMSE values. If \code{box_plot = FALSE}, only the RMSE table is returned.
#' @export
#' @import ggpubr mgcv tibble gss
#' @author Shijia Li
#'
#' @examples
#' \dontrun{
#' # Assuming 'design_data' is obtained from the Design function and 'otu_data' contains the pre-processed OTU data:
#' cv_results <- k.fold.CV(Data_for_Reg = design_data,
#'                         pre_processed_data = otu_data,
#'                         unique_values = 5,
#'                         poly_degree = 3,
#'                         k = 10,
#'                         box_plot = TRUE,
#'                         sig_label = TRUE)
#'
#' # To view the RMSE table:
#' cv_results[['RMSE table']]
#'
#' # To display the boxplots:
#' cv_results[['Boxplot of CV RMSE']]
#' }
#'
k.fold.CV <- function(Data_for_Reg, pre_processed_data, unique_values, poly_degree = 3,k = 10, box_plot = TRUE, sig_label = FALSE) {
  set.seed(123)
  n <- nrow(Data_for_Reg)
  folds <- sample(rep(1:k, length.out = n))

  # 1. Extract independent variable information from output result of the function Design
  inde_vars = setdiff(colnames(Data_for_Reg),colnames(pre_processed_data))
  inde_vars_data = Data_for_Reg[inde_vars]
  rmse_list = list()
  cv_results = list()

  for (i in colnames(inde_vars_data[,-c(1,ncol(inde_vars_data))])){
    Sub_Data <- cbind(pre_processed_data, inde_vars_data[1], inde_vars_data[i])
    Sub_Data <- Sub_Data[Sub_Data[,i] == 1,]
    rmse_df = data.frame(glm_rmse = numeric(),ns_rmse = numeric(),ss_rmse = numeric(),lo_rmse = numeric(),stringsAsFactors = FALSE)

    for (j in colnames(pre_processed_data)){
      # Extract OTU j and the corresponding independent variable
      Sub_data_for_fit = cbind(Sub_Data[j],Sub_Data[ncol(Sub_Data)-1])
      Sub_data_filout = Sub_data_for_fit
      GCV = rep(0,5)
      if (length(unique(Sub_data_filout[,1])) > unique_values){
        rmse_glm_list = c()
        rmse_ns_list = c()
        rmse_ss_list = c()
        rmse_lo_list = c()
        for (z in 1:k) {
          test_indices = which(folds == z)
          train_data = na.omit(Sub_data_filout[-test_indices, ])
          test_data = na.omit(Sub_data_filout[test_indices, ])

          test_data <- test_data[test_data$Time >= min(train_data$Time) &
                                   test_data$Time <= max(train_data$Time), ]

          if (nrow(test_data) == 0) {
            warning("Test data is empty after filtering, skipping this fold.")
            next
          }

          ### glm model
          glm_formula = as.formula(paste(j, "~ poly(Time,degree = ", poly_degree, ")"))
          glm_model = glm(formula = glm_formula, family = gaussian, data = train_data)
          glm_model_step = step(glm_model, direction = 'both', trace = 0)
          glm_predictions = predict(glm_model_step, test_data)

          ### ns model
          for (d in 1:5) {
            ns_formula = as.formula(paste(j, " ~ splines::ns(Time, df = ", 2 + d, ")", sep = ""))
            ns_model = mgcv::gam(formula = ns_formula, data = train_data)
            GCV[d] = summary(ns_model)$sp.criterion[[1]]
          }
          best_knots = c(1:5)[which.min(GCV)]
          knots <- attr(splines::ns(train_data$Time, df = 2 + best_knots), "knots")
          best_formula = as.formula(paste(j, " ~ splines::ns(Time, df = ", 2 + best_knots, ")", sep = ""))
          best_ns_model = mgcv::gam(formula = best_formula, data = train_data)
          ns_predictions = predict(best_ns_model, test_data)

          ### smooth spline model
          ss_formula = as.formula(paste(j, "~ Time"))
          ss_model = gss::ssanova(ss_formula, data = train_data)
          if (!is.null(ss_model)) {
            ss_predictions = predict(ss_model, test_data)
          } else {
            ss_predictions <- rep(NA, nrow(test_data))
          }

          ### loess model
          lo_formula = as.formula(paste(j, "~ Time"))
          lo_model = loess(lo_formula, data = train_data)
          lo_predictions = predict(lo_model, test_data)

          true_values = test_data[[j]]
          rmse_glm = sqrt(mean((glm_predictions - true_values)^2, na.rm = TRUE))
          rmse_ns = sqrt(mean((ns_predictions - true_values)^2, na.rm = TRUE))
          rmse_ss = sqrt(mean((ss_predictions - true_values)^2, na.rm = TRUE))
          rmse_lo = sqrt(mean((lo_predictions - true_values)^2, na.rm = TRUE))

          rmse_glm_list = c(rmse_glm_list, rmse_glm)
          rmse_ns_list = c(rmse_ns_list, rmse_ns)
          rmse_ss_list = c(rmse_ss_list, rmse_ss)
          rmse_lo_list = c(rmse_lo_list, rmse_lo)
        }
      } else {
        next
      }
      mean_rmse_glm = mean(rmse_glm_list)
      mean_rmse_ns = mean(rmse_ns_list)
      mean_rmse_ss = mean(rmse_ss_list)
      mean_rmse_lo = mean(rmse_lo_list)
      mean_rmse = data.frame(mean_rmse_glm,mean_rmse_ns,mean_rmse_ss,mean_rmse_lo)
      rownames(mean_rmse) = j
      colnames(mean_rmse) = c('maSigPro','MicrobTiSDA','MetaDprof','LOWESS')
      rmse_df = rbind(rmse_df,mean_rmse)
    }
    rmse_list[[i]] = rmse_df
  }

  if (box_plot == TRUE) {
    rmse_plot = list()
    for (i in names(rmse_list)) {
      cv_result = rmse_list[[i]] %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Features') %>%
        pivot_longer(cols = -Features, names_to = 'Method', values_to = "RMSE")
      rmse_plot[[i]] = ggplot(data = cv_result, aes(x = Method, y = RMSE, fill = Method))+
        geom_boxplot()+
        theme_minimal()+
        labs(title = "",y = 'RMSE',x = "")+
        theme(axis.title.y = element_text(size = 5),
              axis.text.y = element_text(size = 3),
              axis.text.x = element_text(size = 5,angle = 30),
              legend.position = "none",
              plot.margin = margin(1, 3, 1, 3))
      if (sig_label) {
        rmse_plot[[i]] = rmse_plot[[i]]+
          stat_compare_means(comparisons = list(c("maSigPro", "MicrobTiSDA"),
                                                c("maSigPro", "MetaDprof"),
                                                c("maSigPro", "LOWESS"),
                                                c("MicrobTiSDA", "MetaDprof"),
                                                c("MicrobTiSDA", "LOWESS"),
                                                c("LOWESS", "MetaDprof")),
                             method = 't.test', label = 'p.signif', step.increase = 0.1)
      }
    }
    cv_results = list(rmse_list,rmse_plot)
    names(cv_results) = c('RMSE table','Boxplot of CV RMSE')
    return(cv_results)
  } else {
    cv_results = rmse_list
    names(cv_results) = c('RMSE table')
    return(cv_results)
  }
}
