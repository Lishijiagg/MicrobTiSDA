#' Interpolate Time-Series Data Based on Sample Time
#'
#' @description
#' This function performs interpolation on a data frame or matrix (e.g., OTU/ASV counts or other time-series measurements)
#'     using corresponding metadata time points. For each unique subject (as defined by a subject ID), the function constructs
#'     a full time series between the minimum and maximum time points and applies interpolation (defaulting to cubic interpolation)
#'     to generate data for missing time points. The function returns both the interpolated time-series data and the associated
#'     updated metadata.
#'
#' @details
#' This functin processes the input data and metadata by interating over each unique subject ID defined in \code{Sample_ID}. For
#'     each subject, it subsets and sorts the metadata by \code{Sample_Time} and constructs a complete time series from the minimum
#'     to maximum time values with a step of 1. It then extracts the corresponding data columns and performs interpolation (Using
#'     the specified \code{interp_method}, with \code{cubic} as the default) on each feature across the full time series. Simultaneously,
#'     updated metadata is generated for the interpolated time points, preserving the subject ID and group information as indicated by
#'     \code{Group_var}. The function returns a list containing the interpolated data matrix and the corresponding updated metadata.
#'
#' @param Data A data frame where rows represent OTUs/ASVs and columns represent samples Or the output of the function
#'     \code{\link[MicrobTiSDA]{Data.filter}}.
#' @param metadata A data frame. Containing information about all samples, including at least the grouping of all samples as well as
#'     individual information (\code{Group} and \code{ID}), the sampling \code{Time} point for each sample, and other relevant information.
#' @param Sample_Time A character string specifying the column name in \code{metadata} that contains time information.
#' @param Sample_ID A character string specifying the column name in \code{metadata} that identifies unique samples of each subject.
#' @param interp_method A character string specifying the interpolation method to be used by \code{\link[pracma]{interp1}}.
#'                      Default is \code{'cubic'}. Other methods accepted by \code{interp1} (e.g., \code{'linear'})
#'                      can also be used.
#' @param Group_var A character string specifying the column name in \code{metadata} that indicates group membership.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{Interpolated_Data}{A data frame containing the interpolated data. Columns represent interpolated time points
#'                            (labeled with the sample ID and time), and rows correspond to features.}
#'   \item{Interpolated_Data_metadata}{A data frame containing metadata for the interpolated data, including the time points,
#'                                       sample IDs, and group information.}
#' }
#'
#' @importFrom pracma interp1
#' @author Shijia Li
#' @export
#'
#' @examples
#' \dontrun{
#' # Example data: 5 features across 8 samples with time points from two subjects.
#' set.seed(123)
#' Data <- matrix(sample(1:100, 40, replace = TRUE), nrow = 5)
#' rownames(Data) <- paste0("Feature", 1:5)
#' colnames(Data) <- paste0("Sample", 1:8)
#'
#' # Create metadata with time points, sample IDs, and group assignments.
#' metadata <- data.frame(
#'   Time = c(1, 3, 5, 7, 2, 4, 6, 8),
#'   ID = c(rep("Subject1", 4), rep("Subject2", 4)),
#'   Group = c(rep("A", 4), rep("B", 4)),
#'   row.names = paste0("Sample", 1:8)
#' )
#'
#' # Interpolate the data using cubic interpolation.
#' interp_results <- Data.interpolate(Data = Data,
#'                                    metadata = metadata,
#'                                    Sample_Time = "Time",
#'                                    Sample_ID = "ID",
#'                                    interp_method = "cubic",
#'                                    Group_var = "Group")
#' }
Data.interpolate = function(Data,
                            metadata,
                            Sample_Time,
                            Sample_ID,
                            interp_method = 'cubic',
                            Group_var){

  sample_id = unique(metadata[[Sample_ID]])
  interpolated_data = data.frame(matrix(ncol = 0,nrow = nrow(Data)))
  interpolated_meta = data.frame()

  for (i in sample_id) {
    meta_id = subset(metadata,metadata[[Sample_ID]] == i)
    sorted_meta_id = meta_id[order(meta_id[[Sample_Time]]),]
    Group = unique(sorted_meta_id[[Group_var]])
    full_time_series = seq(min(sorted_meta_id[,Sample_Time]),max(sorted_meta_id[,Sample_Time]),1)
    data_id = Data[,rownames(sorted_meta_id)]

    intp_data = as.data.frame(matrix(0,nrow = nrow(data_id),ncol = length(full_time_series)))
    rownames(intp_data) = rownames(data_id)
    colnames(intp_data) = paste(i,"_T_",full_time_series,sep = "")

    intp_meta = as.data.frame(matrix(0,nrow = ncol(intp_data),ncol = 0))
    rownames(intp_meta) = colnames(intp_data)
    intp_meta$Time = full_time_series
    intp_meta$ID = rep(i,nrow(intp_meta))
    intp_meta$Group = rep(Group,nrow(intp_meta))

    for (j in rownames(data_id)) {
      intp_data[j,] = pracma::interp1(x = as.numeric(sorted_meta_id[,Sample_Time]),
                                      y = as.numeric(data_id[j,]),
                                      xi = full_time_series,
                                      method = interp_method)
    }

    interpolated_data = cbind(interpolated_data,intp_data)
    interpolated_meta = rbind(interpolated_meta,intp_meta)
  }

  interp_results = list(interpolated_data,interpolated_meta)
  names(interp_results) = c('Interpolated_Data','Interpolated_Data_metadata')
  return(interp_results)
}
