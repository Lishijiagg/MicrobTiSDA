
#' @title Fit Mixed-Effects Spline Regression Models with GAM for Microbial Features on Group-Level
#' @description
#' This function fits mixed-effect spline regression models using GAM to capture nonlinear temporal trends in microbial
#'     community data at the group level. It incorporates random effects for individual IDs within each group and uses
#'     natural splines to model the relationship between time and microbial feature abundances.
#' @details
#' This function fits mixed-effects spline regression models. Unlike \code{\link[MicrobTiSDA]{Reg.SPLR}}, this function
#'     captures the overall temporal dynamics of microbial features across different groups. Assuming the dataset contains
#'     \eqn{\textit{m}} groups, each with \eqn{\textit{n}} subjects, the model is formulated as:
#'     \deqn{x_{(mi)}(t) = \beta_{m0} + \sum_{k=1}^{K} \beta_{mk} \cdot N_{mk}(t) + b_{n(m)} + \epsilon}
#'     where \eqn{x_{(mi)}(t)} represents the abundance of microbial feature \eqn{\textit{i}} at time point \eqn{\textit{t}}
#'     in group \eqn{\textit{m}}. The random effect \eqn{b_{n(m)}} reflects the departure of individual \eqn{\textit{n}} in group
#'     \eqn{\textit{m}} from overall population average effects. The parameter \eqn{K} refers to the number of basis functions,
#'     which is equal to the number of knots plus one. (i.e., \eqn{K = number of knots + 1}).
#'
#' The \code{Reg.MESR} function first extracts independent variables from a design matrix (produced by the \code{\link[MicrobTiSDA]{Design}})
#'     by removing columns corresponding to the pre-processed OTU data. For each design dummy variable (excluding the first
#'     and last columns), the function subsets the data to include only the observations where the dummy variable equals 1.
#'     Then, for each OTU in the pre-processed data, it optionally filters out outliers based on a specified z-score threshold.
#'     If the number of unique transformed OTU values exceeds a given threshold (\code{unique_values}), a mixed-effects spline
#'     regression model is fitted using a natural spline on the time variable along with a random effect for the sample ID.
#'     When the \code{Knots} parameter is not provided, the function iterates over a range of knot numbers (from 1
#'     to \code{max_Knots}), selects the optimal model by minimizing the Generalized Cross-Validation (GCV) criterion,
#'     and extracts the corresponding knot locations. Alternatively, if \code{Knots} is provided, it is directly used in
#'     model fitting. The resulting fitted models and associated knots information are organized into nested lists and returned.
#'
#'
#' @param Data_for_Reg The data frame output from the \code{\link[MicrobTiSDA]{Design}}.
#' @param pre_processed_data The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
#'     pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.
#' @param unique_values An integer specifying the minimum number of unique values required in an OTU for spline fitting (default: 5).
#' @param z_score A numeric value specifying the z-score threshold for outlier filtering; if \code{NA}, no
#'     outlier filtering is performed (default: \code{NA}).
#' @param Knots An optional numeric vector specifying the knots to use in the spline regression. If \code{NULL}, the optimal number of
#'     knots is determined by minimizing the GCV criterion (default: \code{NULL}).
#' @param max_Knots An integer indicating the maximum number of knots to consider when selecting the optimal spline model (default: 3).
#' @import mgcv splines
#' @return A list with two elements: \code{fitted_model} is a nested list containing the fitted mixed-effects GAM models for each design
#'     dummy variable and OTU; \code{knots_info_each_model} is a corresponding nested list with the knots used for each model.
#' @export
#' @author Shijia Li
#' @examples
#' \dontrun{
#' # Assuming design_data is obtained from the Design function, otu_data contains pre-processed OTU data,
#' # and metadata is available:
#' fit_result <- Reg.MESR(Data_for_Reg = design_data,
#'                        pre_processed_data = otu_data,
#'                        metadata = metadata,
#'                        unique_values = 5,
#'                        z_score = 2,
#'                        Knots = NULL,
#'                        max_Knots = 5)
#' # Access the fitted model for a specific group and OTU:
#' model_example <- fit_result$fitted_model[['Group1']][['OTU1']]
#' }
Reg.MESR = function(Data_for_Reg,
                    pre_processed_data,
                    unique_values = 5,
                    z_score = NA,
                    Knots = NULL,
                    max_Knots = 3) {

  # 1. Extract independent variables from the output result of the Design function
  inde_var = setdiff(colnames(Data_for_Reg),colnames(pre_processed_data))
  inde_var_data = Data_for_Reg[inde_var]

  fit_data = list()
  knots_info = list()

  # 2. Fit polynomial regression model for each OTU and store the results to the fit_data list
  for (i in colnames(inde_var_data[,-c(1,ncol(inde_var_data))])) {

    set.seed(123)
    GCV = rep(0,max_Knots)
    Sub_data = cbind(pre_processed_data,inde_var_data[1],inde_var_data[i],inde_var_data[ncol(inde_var_data)])
    Sub_data = Sub_data[Sub_data[,i] == 1,]

    for (j in colnames(pre_processed_data)) {

      # Extract OTU j and the corresponding independent variables
      Sub_data_for_fit = cbind(Sub_data[j],Sub_data[ncol(Sub_data)-2],Sub_data[ncol(Sub_data)])

      # Apply z-score filtering
      if (!is.na(z_score)) {

        # calculate Z-score for each OTU
        z_scores = scale(Sub_data_for_fit[1])
        outliers_idx = which(abs(z_scores) > z_score)
        Sub_data_for_fit = Sub_data_for_fit[-outliers_idx,]
      } else {
        Sub_data_for_fit = Sub_data_for_fit
      }

      Sub_data_for_fit$ID = as.factor(Sub_data_for_fit$ID)

      # Fit mixed-effect Spline regression models
      if (length(unique(Sub_data_for_fit[,1])) > unique_values) {

        if (is.null(Knots)) {
          for (k in 1:max_Knots) {
            formula = as.formula(paste(j," ~ splines::ns(Time, df = ",1+k,")+s(ID,bs = 're')",
                                       sep = ""))
            sp_model = mgcv::gam(formula = formula,
                                  data = Sub_data_for_fit,
                                  method = "GCV.Cp")
            GCV[k] = sp_model$gcv.ubre
          }
          best_knots = c(1:max_Knots)[which.min(GCV)]
          knots = attr(splines::ns(Sub_data_for_fit$Time,df = 1+best_knots),"knots")
          knots_info[[i]][[j]] = knots
          best_formula = as.formula(paste(j," ~ splines::ns(Time, df = ",1+best_knots,")+s(ID,bs = 're')",
                                          sep = ""))
          fit_data[[i]][[j]] = mgcv::gam(formula = best_formula,
                                         data = Sub_data_for_fit,
                                         method = 'GCV.Cp')
        } else {
          Knots_str = paste(Knots,collapse = ",")
          formula = as.formula(paste(j," ~ splines::ns(Time, knots = c(",Knots_str,"))+s(ID,bs='re')",
                                     sep = ''))
          fit_data[[i]][[j]] = mgcv::gam(formula = formula,data = Sub_data_for_fit,
                                         method = 'GCV.Cp')
          knots_info[[i]][[j]] = Knots
        }
      } else {
        fit_data[[i]][[j]] = NULL
        knots_info[[i]][[j]] = NULL
        cat("The model of",j,"in group",i,"was excluded due to less than specified unique transformed value for model fitting\n")
      }
    }
  }
  fit_result = list(fit_data,knots_info)
  names(fit_result) = c("fitted_model","knots_info_each_model")
  return(fit_result)
}
