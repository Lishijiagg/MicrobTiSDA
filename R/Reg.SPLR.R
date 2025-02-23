#' @title Fit Spline Regression Models with Knot Selection
#' @description
#' The \code{Reg.SPLR} function fits natural spline regression models to time-series OTU data for each subgroup
#'     defined by design dummy variables. It uses \code{\link[mgcv]{gam}} to model the relationship between OTU
#'     abundance and time, incorporating a z-score based outlier filtering step (optional) and selecting the optimal number
#'     of knots via the Generalized Cross-Validation (GCV) criterion when not provided by the user.
#' @details
#' \code{Reg_SPLR} utilize Generalised Additive Model (GAM,see \code{\link[mgcv]{gam}}) to fit natural splines. Here, let \eqn{x_i(t)} denote
#'     the transformed abundance of microbial feature \eqn{i} at time point \eqn{t}. To explain the evolution of microbial abundance along the
#'     time (\eqn{t}), the following generalized additive model is considered
#'     \deqn{x_i(t) = \beta_0 + \sum_{k=1}^K \beta_k \cdot N_k(t) + \epsilon}
#'     where \eqn{\beta_0} is the intercept term, \eqn{\beta_k} denotes the regression coefficients for the k-th natural spline basis function
#'     \eqn{N_k(t)}, and \eqn{\epsilon} is the error term. The parameter \eqn{K} refers to the number of basis functions, which is equal to the
#'     number of knots plus one. (i.e., \eqn{K = number of knots + 1}).
#'
#' If a z-score threshold is specified via \code{z_score}, observations with absolute z-scores exceeding this threshold are removed prior to model
#'     fitting. When the OTU has more unique values than specified by \code{unique_values}, the function searches over a range of knot numbers
#'     (from 1 up to \code{max_Knots}) to select the optimal model based on the GCV criterion; if \code{Knots} is provided, it is used directly
#'     in constructing the natural spline basis. The fitted models and their corresponding knots information are stored in a nested list
#'     structure and returned.
#'
#'
#' @param Data_for_Reg The data frame output from the \code{\link[MicrobTiSDA]{Design}}.
#' @param pre_processed_data The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
#'     pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.
#' @param z_score A numeric value specifying the z-score threshold for outlier filtering; if \code{NA}, no
#'     outlier filtering is performed (default: \code{NA}).
#' @param unique_values An integer specifying the minimum number of unique values required in an OTU for spline fitting (default: 5).
#' @param Knots An optional numeric vector specifying the knots to use in the spline regression. If \code{NULL}, the optimal number of
#'     knots is determined by minimizing the GCV criterion (default: \code{NULL}).
#' @param max_Knots An integer indicating the maximum number of knots to consider when selecting the optimal spline model (default: 5).
#'
#' @return A list with two elements: \code{fitted_model} is a nested list of the fitted spline regression models for each design dummy
#'     variable and OTU; \code{knots_info_each_model} contains the corresponding knots information for each model.
#' @export
#' @import mgcv splines
#'
#' @examples
#' \dontrun{
#' # Assuming Data_for_Reg is obtained from a design step and pre_processed_data contains OTU data:
#' result <- Reg.SPLR(Data_for_Reg, pre_processed_data, z_score = 2, unique_values = 5, Knots = NULL, max_Knots = 5)
#' # Access the fitted model for a particular design dummy variable and OTU:
#' fitted_model_example <- result$fitted_model[['Group1']][['OTU_name']]
#' }
Reg.SPLR = function(Data_for_Reg,
                    pre_processed_data,
                    z_score = NA,
                    unique_values = 5,
                    Knots = NULL,
                    max_Knots = 5) {

  # 1. Extract independent variable information from output result of design function
  inde_vars = setdiff(colnames(Data_for_Reg),colnames(pre_processed_data))
  inde_vars_data = Data_for_Reg[inde_vars]
  fit_data = list()
  knots_info = list()

  # 2. Fit spline regression model of each OTU and store the results to a list named fit_data
  for (i in colnames(inde_vars_data[,-c(1,ncol(inde_vars_data))])) { # iterate through design dummy variables
    set.seed(123)
    GCV = rep(0,max_Knots)
    Sub_Data <- cbind(pre_processed_data, inde_vars_data[1], inde_vars_data[i])
    Sub_Data <- Sub_Data[Sub_Data[,i] == 1,]

    for (j in colnames(pre_processed_data)) {
      Sub_Data_for_fit = cbind(Sub_Data[j], Sub_Data[ncol(Sub_Data)-1])

      # Apply Z-score filtering
      if (!is.na(z_score)) {
        # calculate Z-score for each OTU
        z_scores = scale(Sub_Data_for_fit[1])
        outliers_idx = which(abs(z_scores) > z_score)
        Sub_data_filout = Sub_Data_for_fit[-outliers_idx,]
      } else {
        Sub_data_filout = Sub_Data_for_fit
      }

      # fit natural spline regression for each filtered OTU
      if (length(unique(Sub_data_filout[,j])) > unique_values) {

        if (is.null(Knots)) {
          for (k in 1:max_Knots) {
            formula = as.formula(paste(j," ~ splines::ns(Time, df = ",1+k,")",sep = ""))
            #formula = as.formula(paste(j," ~ ns(Time, knots = ",k,")",sep = ""))
            sp_model = mgcv::gam(formula = formula,data = Sub_data_filout)
            GCV[k] = summary(sp_model)$sp.criterion[[1]]
          }
          best_knots = c(1:max_Knots)[which.min(GCV)]
          knots <- attr(splines::ns(Sub_data_filout$Time, df = 1+best_knots), "knots")
          knots_info[[i]][[j]] = knots
          best_formula = as.formula(paste(j, " ~ splines::ns(Time, df = ",1+best_knots,")",
                                          sep = ""))
          fit_data[[i]][[j]] = mgcv::gam(formula = best_formula,data = Sub_data_filout)
        } else {
          Knots_str = paste(Knots,collapse = ",")
          formula = as.formula(paste(j," ~ splines::ns(Time, knots = c(",Knots_str,"))",sep = ''))
          fit_data[[i]][[j]] = mgcv::gam(formula = formula,data = Sub_data_filout)
          knots_info[[i]][[j]] = Knots
        }

      } else {
        fit_data[[i]][[j]] = NULL
        knots_info[[i]][[j]] = NULL
      }
    }
  }

  fit_result = list(fit_data,knots_info)
  names(fit_result) = c("fitted_model","knots_info_each_model")

  # 3. Return fitted models
  return(fit_result)
}
