#' @title Transform Microbial Composition Data.
#'
#' @description
#' This function applies the modified centered log-ratio (MCLR) transformation function \code{\link[MicrobTiSDA]{mclr.transform}}
#'     to a data matrix (e.g., OTU/ASV counts). When a grouping variable is provided the transformation is applied separately for
#'     each group defined in the metadata.
#' @details
#' The function transforms the input data using the MCLR method. If no grouping variable is provided (i.e. \code{Group_var}
#'     is \code{NULL}), the transformation is applied to the entire dataset. If a single grouping variable is specified,
#'     the data is partitioned into subsets corresponding to the unique groups in the metadata, and the transformation
#'     is applied to each subset separately; the results are then combined using row binding. For multiple grouping variables,
#'     a composite grouping factor is created using the interaction of the specified variables, and the transformation is
#'     applied to each composite group in a similar manner.
#'
#' @param Data A data frame or matrix of microbial compositional data, with rows representing microbial features (OTUs/ASVs)
#'     and columns representing samples.
#' @param metadata A data frame. Containing information about all samples, including at least the grouping of all samples as well as
#'     individual information (\code{Group} and \code{ID}), the sampling \code{Time} point for each sample, and other relevant information.
#' @param Group_var A string or a vector. This specifies the grouping variables, which should match the column names in
#'     the \code{metadata} used to designate sample groups, and for pre-processing OTU data of each group or individual separately.
#'     For instance, to split the OTU table based on the \code{"Group"} variable, set \code{Group_var = "Group"};
#'     to split the data based on the \code{"Group"} and \code{"Diet"} (if in \code{metadata})categorical variables to study the
#'     interaction between different grouping variables, set \code{Group_var = c("Group","Diet")}.
#'
#' @return A data frame containing the MCLR-transformed data with the same structure as the input microbial compositional data.
#' @author Shijia Li
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data matrix (5 features x 10 samples)
#' set.seed(123)
#' Data <- matrix(sample(1:100, 50, replace = TRUE), nrow = 5)
#' rownames(Data) <- paste0("Feature", 1:5)
#' colnames(Data) <- paste0("Sample", 1:10)
#'
#' # Create example metadata with a grouping variable
#' metadata <- data.frame(Group = rep(c("A", "B"), each = 5))
#' rownames(metadata) <- paste0("Sample", 1:10)
#'
#' # Apply MCLR transformation to the entire dataset
#' transformed_data <- Data.trans(Data, metadata, Group_var = NULL)
#'
#' # Apply MCLR transformation separately for each group
#' transformed_data_by_group <- Data.trans(Data, metadata, Group_var = "Group")
#' }
#'
Data.trans = function(Data,metadata,Group_var) {

  if (is.null(Group_var)) {

    Transformed_data = as.data.frame(t(mclr.transform(Data)))

  } else if (length(Group_var) < 2) {

    groups = unique(metadata[,Group_var])
    group_OTU = list()
    for (group in groups) {
      samples_in_group = rownames(subset(metadata, metadata[,Group_var] == group))
      group_data = Data[,samples_in_group]
      Transformed_Data = as.data.frame(t(mclr.transform(group_data)))
      group_OTU[[group]] = Transformed_Data
    }
    transformed_data = dplyr::bind_rows(group_OTU)

  } else {

    groups_combo = metadata[, Group_var, drop = FALSE]
    metadata$Group_combo = interaction(metadata[,colnames(groups_combo)])

    groups = unique(metadata[,"Group_combo"])
    group_OTU = list()
    for (group in groups) {
      samples_in_group = rownames(subset(metadata, metadata[,Group_combo] == group))
      group_data = Data[,samples_in_group]

      Transformed_Data = as.data.frame(t(mclr.transform(group_data)))
      group_OTU[[group]] = Transformed_Data
    }
    transformed_data = dplyr::bind_rows(group_OTU)
  }

  return(transformed_data)
}
