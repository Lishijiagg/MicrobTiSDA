#' @title Cluster OTU Time-Series Data Based on Regression Model prediction and Generate Dendrogram Plots
#'
#' @description
#' This function performs hierarchical clustering on predicted OTU time-series data for different groups
#'     and generates corresponding dendrogram plots. For each group in the input list, the function computes a
#'     correlation-based distance matrix, performs hierarchical clustering using the specified clustering method
#'     (e.g. \code{average}), and then converts the result into a dendrogram.
#'
#' @details
#' For each group in the input \code{predicted_data}, the function first extracts the predicted OTU data (excluding the
#'     last column, which is assumed to contain time information) and computes a correlation matrix, which is converted
#'     into a distance matrix via
#'     \deqn{d_{\text{corr}}(x,y) = 1-\frac{{\sum_{i=1}^{n}(x_i-\bar{y})}}{{\sqrt{{\sum_{i=1}^{n}(x_i-\bar{x})^2}} \sqrt{{\sum_{i=1}^{n}(y_i-\bar{y})^2}}}}}
#'     where \eqn{x} and \eqn{y} represent the two OTU time series being compared, \eqn{n} denotes the total number of time points, and
#'     \eqn{\bar{x}} and \eqn{\bar{y}} denote the means of the respective time series. Hierarchical clustering is
#'     performed on the above distance matrix using the method specified in \code{clust_method}.
#'
#'
#'
#' @param predicted_data The output data frame from the \code{\link[MicrobTiSDA]{Pred.data}}.
#' @param clust_method A string, the agglomeration method to be used. This argument should be one of "ward.D", "ward.D2", "single",
#'     "complete", "average", "mcquitty", "median", "centroid". Detail see \code{\link[stats]{hclust}}.
#' @param font_size A numeric value specifying the font size for text labels in the dendrogram plots (default: \code{0.2}).
#' @param dend_title_size A numeric value specifying the font size of the dendrogram plot title (default: \code{15}).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{predicted_data}{The original input list of predicted data.}
#'   \item{cluster_results}{A list of hierarchical clustering objects (one per group).}
#'   \item{cluster_figures}{A list of ggplot2 objects containing the dendrogram plots for each group.}
#' }
#'
#' @export
#' @author Shijia Li
#' @importFrom ggplot2 ggplot
#' @importFrom cluster silhouette
#' @importFrom ggdendro dendro_data
#' @importFrom dendextend color_branches
#'
#' @examples
#' \dontrun{
#' # Assuming you have a list of predicted data for each group (my_predicted_data) generated by \code{\link[MicrobTiSDA]{Pred.data}}
#' # and you wish to use the "average" linkage method:
#' result <- Data.cluster(predicted_data = my_predicted_data,
#'                        clust_method = "average",
#'                        font_size = 0.2,
#'                        dend_title_size = 15)
#'
#' # To view the dendrogram plot for a particular group:
#' print(result$cluster_figures[["Group1"]])
#' }
Data.cluster = function(predicted_data,clust_method,font_size=0.2,dend_title_size=15) {

  cluster_results = list()
  cluster_figures = list()

  for (i in names(predicted_data)) {

    pred_data = as.data.frame(predicted_data[[i]][, -ncol(predicted_data[[i]])])

    if (is.null(pred_data) || (!is.data.frame(pred_data) && !is.matrix(pred_data)) || ncol(pred_data) < 3) {
      cluster_results[[i]] = NULL
      cat("Group", i, "was excluded due to having fewer than 3 OTUs, making it unsuitable for clustering.\n")
    } else {

      cor_matrix = cor(pred_data)
      pred_data_dist = as.dist((1 - cor_matrix))

      pred_data_dist_hc = hclust(pred_data_dist, clust_method)

      dend = as.dendrogram(pred_data_dist_hc)
      cluster_figures[[i]] = ggplot() +
        geom_segment(data = ggdendro::dendro_data(dend, type = "rectangle")$segments,
                     aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(data = ggdendro::dendro_data(dend, type = "rectangle")$labels,
                  aes(x = x, y = y, label = label),
                  size = font_size, hjust = 0) +
        coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
        labs(title = paste("OTU profiles clustering plot for group", i)) +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_blank(),
              plot.title = element_text(size = dend_title_size))
      print(cluster_figures[[i]])
      cluster_results[[i]] = pred_data_dist_hc
    }
  }
  cluster_data = list(predicted_data,cluster_results, cluster_figures)
  names(cluster_data) = c("predicted_data","cluster_results", "cluster_figures")
  return(cluster_data)
}

