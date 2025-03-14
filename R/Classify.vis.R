#' Visualize Microbial Feature Classification Results Using NMDS
#'
#' @description
#' This function generates a non-metric multidimensional scaling (NMDS) plot to visualize OTU/ASV classification
#'     results. It combines an OTU/ASV table with predicted group labels to produce a scatter plot in NMDS space,
#'     where each sample is colored according to its predicted group.
#'
#' @details
#' The function expects as input the results of output of the function \code{\link[MicrobTiSDA]{Data.rf.classifier}}.
#'     It conputes NMDS coordinates based on the OTU/ASV data using the specified distance method (defaulting to
#'     Bray-Curtis) and then maps the predicted group labels onto the NMDS coordinates. The visualization is produced
#'     using \code{ggplot2}, with polygons delineating the convex hull of each predicted group.
#'
#'
#' @param classified_results A list object of \code{\link[MicrobTiSDA]{Data.rf.classifier}} output.
#' @param dist_method Dissimilarity index, defaults to \code{bray}. For other options, including
#'     \code{manhattan}, \code{euclidean}, \code{canberra} etc., see \code{\link[vegan]{vegdist}}.
#' @param fig_title A character string specifying the title of the NMDS plot. Default is \code{'NMDS visualization'}.
#' @param legd_title A character string for the legend title representing the predicted groups. Default is \code{'Predicted Groups'}.
#' @param points_size A numeric value specifying the size of the points in the plot. Default is \code{1}.
#' @param legend_title_size A numeric value specifying the font size of the legend title. Default is \code{8}.
#' @param legend_text_size A numeric value specifying the font size of the legend text. Default is \code{6}.
#' @param axis_title_size A numeric value specifying the font size of the axis titles. Default is \code{6}.
#' @param axis_text_size A numeric value specifying the font size of the axis text. Default is \code{6}.
#'
#' @return A \code{ggplot2} object displaying the NMDS plot with samples colored by their predicted group and convex hulls
#'         outlining each group.
#'
#' @importFrom vegan metaMDS
#' @importFrom plyr ddply
#' @importFrom ggplot2 ggplot
#' @author Shijia Li
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'classified_results' is the output from a classification function:
#' nm_plot <- Classify.vis(classified_results = results,
#'                         dist_method = 'bray',
#'                         fig_title = 'NMDS Plot',
#'                         legd_title = 'Predicted Groups',
#'                         points_size = 1.5)
#' }
#'
Classify.vis = function(classified_results,
                           dist_method = 'bray',
                           fig_title = 'NMDS visualization',
                           legd_title = 'Predicted Groups',
                           points_size = 1,
                           legend_title_size = 8,
                           legend_text_size = 6,
                           axis_title_size = 6,
                           axis_text_size = 6) {

  otu = classified_results[[1]]
  #otu = otu[which(rowSums(otu) >= min_total_count),]
  nmds = vegan::metaMDS(t(otu),distance = dist_method)
  result = nmds$points
  result = as.data.frame(cbind(result, rownames(result)))

  # obtain the classification result from the output results of the OTU_classifier function
  predict_group = c(classified_results[[2]],classified_results[[3]])
  predict_group = as.character(predict_group[rownames(result)])

  # figure visualization
  colnames(result)[1:3] = c('NMDS1','NMDS2','samples')
  result$NMDS1 = as.numeric(as.character(result$NMDS1))
  result$NMDS2 = as.numeric(as.character(result$NMDS2))
  result$samples = as.character(result$samples)
  result = cbind(result, predict_group)

  p = ggplot2::ggplot(result,aes(NMDS1,NMDS2,color = predict_group)) +
    geom_polygon(data = plyr::ddply(result,'predict_group',
                                    function(df) df[chull(df[[1]], df[[2]]),]),fill = NA) +
    labs(title = fig_title,
         color = legd_title) +
    theme(legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size))+
    geom_point(size = points_size)+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'))

  print(p)
  return(p)
}
