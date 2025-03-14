#' @title Select Biomarkers Based on Random Forest Cross-Validation Results
#'
#' @description
#' This function extracts the top biomarkers from a random forest classification result based on cross-validation
#'     and user specified number of microbial features. It updates the cross-validation plot by adding a vertical
#'     dashed line at the specified number of features, and then selects the top features (biomarkers) based on
#'     their importance ranking.
#'
#' @details
#' The function takes an object (usually the output from \code{\link[MicrobTiSDA]{Data.rf.classifier}}, which includes
#'     a cross-validation plot, an OTU importance table, and the original input data) and a user-specified
#'     number of features to select. It then updates the cross-validation plot by adding a vertical dashed line
#'     at the position corresponding to the number of selected features. Next, it extracts the top features from the
#'     OTU importance table (ordered by Mean Decrease Accuracy) and creates q table of these features
#'     from the original microbial feature table. The function returns a list that includes both the transposed biomarker
#'     table and the modified cross-validation plot.
#'
#'
#' @param rf A list containing the results of the random forest classification. Default to \code{\link[MicrobTiSDA]{Data.rf.classifier}}.
#' @param feature_select_num A numeric value specifying the number of top features (biomarkers) to select. Typically, the numer of
#'     specified biomarkers needs to be determined by the user based on the cross-validation result plot output by
#'     \code{\link[MicrobTiSDA]{Data.rf.classifier}}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{OTU_importance}{A data frame of the selected biomarkers (transposed feature table).}
#'     \item{cross_validation_fig}{A ggplot object of the cross-validation plot with a vertical dashed line indicating the feature selection cutoff.}
#'   }
#' @export
#' @importFrom ggplot2 ggplot
#' @author Shijia Li
#'
#' @examples
#' \dontrun{
#' # Assuming rf_results contains the necessary components from a random forest classification
#' # and you wish to select the top 20 features:
#' result <- Rf.biomarkers(rf = rf_results, feature_select_num = 20)
#' # View the biomarker table
#' print(result$OTU_importance)
#' # View the updated cross-validation plot
#' print(result$cross_validation_fig)
#' }
Rf.biomarkers = function(rf = rf_results,feature_select_num) {
  p = rf$cross_validation +
    geom_vline(xintercept = as.numeric(feature_select_num),linetype = 'dashed')+
    scale_x_continuous(breaks = as.numeric(feature_select_num))

  train_predict = rf$Predicted_results_on_train_set
  test_predict = rf$Predicted_results_on_test_set

  print(p)

  selected_biomarker = rownames(rf$OTU_importance)[1:as.numeric(feature_select_num)]
  feature_table = rf$Input_data
  important_otu_table = feature_table[selected_biomarker]

  output_biomarker_table = as.data.frame(t(important_otu_table))

  output_results = list(output_biomarker_table,train_predict,test_predict,p)
  names(output_results) = c('OTU_importance',
                            'Predicted_results_on_train_set',
                            'Predicted_results_on_test_set',
                            'cross_validation_fig')

  return(output_results)
}
