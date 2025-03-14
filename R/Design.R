#' @title Create Design Matrix for Regression Analysis
#' @description
#' \code{Design} creates the design matrix of dummies for fitting regression models of microbiota in time.
#'
#' @details
#' The main functionality of \code{Design} is to add user-selected sample information to the pre-processed OTU/ASV table
#'     as independent variables for fitting the OTU time series regression models. One necessary independent variable
#'     for fitting is \code{Time}, so the default output of this function is the transformed OTU/ASV table with added
#'     sample Time information. If the user also inputs other qualitative variables such grouping, gender, etc., the
#'     function will define dummy variables to distinguish each group based on the number of qualitative variables
#'     entered by the user and the grouping situation of samples based on qualitative variables. Moreover, the subject ID
#'     of each sample will be added as a column to the generated design matrix.
#'
#' @param metadata A data frame containing information for all samples, which should be identical to the \code{metadata}
#'     received by other functions in \code{MicrobTiSDA}.
#' @param Group_var A string or a vector. Same as the \code{Group_var} in \code{\link[MicrobTiSDA]{Data.trans}}.
#' @param Pre_processed_Data The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
#'     pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.
#' @param Sample_Time A character string indicating the column name in \code{metadata} that contains sample time information.
#' @param Sample_ID A character string indicating the column name in \code{metadata} that contains sample identifiers.
#'
#' @return A data frame that combines the pre-processed data, the standardized time variable (\code{'Time'}), and,
#'     if grouping variables are provided, the design matrix representing unique group interactions, along with the
#'     sample identifier (\code{'ID'}).
#' @export
#' @author Shijia Li
#'
#' @examples
#' \dontrun{
#' # Example metadata with grouping variables
#' metadata <- data.frame(
#'   TimePoint = c(1, 2, 3, 4),
#'   Sample = c('S1', 'S2', 'S3', 'S4'),
#'   GroupA = c('A', 'A', 'B', 'B'),
#'   GroupB = c('X', 'Y', 'X', 'Y')
#' )
#'
#' # Example pre-processed data (e.g., transformed abundance data)
#' Pre_processed_Data <- data.frame(
#'   Feature1 = rnorm(4),
#'   Feature2 = rnorm(4)
#' )
#'
#' # Create design matrix using grouping variables
#' design_data <- Design(metadata, Group_var = c('GroupA', 'GroupB'), Pre_processed_Data,
#'                       Sample_Time = 'TimePoint', Sample_ID = 'Sample')
#'
#' # Create design data without grouping variables
#' design_data_no_group <- Design(metadata, Group_var = NULL, Pre_processed_Data,
#'                                Sample_Time = 'TimePoint', Sample_ID = 'Sample')
#' }
Design <- function(metadata, Group_var = NULL, Pre_processed_Data, Sample_Time, Sample_ID) {

  colnames(metadata)[colnames(metadata) == Sample_Time] = 'Time'
  colnames(metadata)[colnames(metadata) == Sample_ID] = 'ID'

  # Extract sample time information from metadata
  Time = metadata$Time
  ID = metadata$ID

  if (!is.null(Group_var)) {
    # Extract columns corresponding to the categorical variables in metadata
    selected_columns <- metadata[, Group_var, drop = FALSE]

    # Convert the extracted column of categorical variables into factors
    selected_columns <- lapply(selected_columns, factor)

    # Create interaction terms to represent combinations of grouping variables
    interaction_terms <- interaction(selected_columns, drop = TRUE)

    # Create the design matrix
    design_matrix <- model.matrix(~interaction_terms - 1)

    # Name each column
    colnames(design_matrix) <- levels(interaction_terms)

    # combine time, design matrix, and transformed data
    Data_for_Reg <- as.data.frame(design_matrix) %>%
      cbind(Time, .) %>%
      cbind(Pre_processed_Data, .)
    Data_for_Reg$ID = ID
  } else {
    Data_for_Reg <- cbind(Pre_processed_Data, Time)
  }

  # return combined
  return(Data_for_Reg)
}
