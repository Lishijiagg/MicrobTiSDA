#' @title Perform Hierarchical Clustering on Predicted Microbial Abundance Data and Visualize Dendrograms
#' @description
#' The \code{Data.cluster} function performs hierarchical clustering on predicted OTU time-series data to group
#'     OTUs with similar temporal patterns. It computes a correlation-based distance matrix and applies hierarchical
#'     clustering using a specified method. Depending on the \code{auto_cutree} parameter, the function either automatically
#'     determines the optimal number of clusters using silhouette analysis or prompts the user to manually specify a
#'     cut-off line. The resulting dendrograms are then annotated and returned along with the clustering objects.
#'
#' @details
#' For each group in the input \code{predicted_data}, the function first extracts the predicted OTU data (excluding the
#'     last column, which is assumed to contain time information) and computes a correlation matrix, which is converted
#'     into a distance matrix via
#'     \deqn{d_{\text{corr}}(x,y) = 1-\frac{{\sum_{i=1}^{n}(x_i-\bar{y})}}{{\sqrt{{\sum_{i=1}^{n}(x_i-\bar{x})^2}} \sqrt{{\sum_{i=1}^{n}(y_i-\bar{y})^2}}}}}
#'     where \eqn{x} and \eqn{y} represent the two OTU time series being compared, \eqn{n} denotes the total number of time points, and
#'     \eqn{\bar{x}} and \eqn{\bar{y}} denote the means of the respective time series.
#'
#' Hierarchical clustering is performed on the above distance matrix using the method specified in \code{clust_method}. When \code{auto_cutree} is \code{TRUE},
#'     the function iteratively cuts the dendrogram for cluster counts ranging from 2 to n-1 (where n is the
#'     number of OTUs) and calculates the average silhouette width for each case, selecting the optimal number of
#'     clusters that maximizes the silhouette width. The dendrogram is then annotated with the suggested optimal
#'     cluster number. If \code{auto_cutree} is \code{FALSE}, the dendrogram is displayed and the user is prompted
#'     to input a cut-off height, which is used to annotate the dendrogram with a dashed line and corresponding label.
#'     In both scenarios, the function returns a list containing the hierarchical clustering results and the
#'     corresponding ggplot2 dendrogram figures.
#'
#' @param predicted_data The output data frame from the \code{\link[MicrobTiSDA]{Pred.data}}.
#' @param clust_method A string, the agglomeration method to be used. This argument should be one of "ward.D", "ward.D2", "single",
#'     "complete", "average", "mcquitty", "median", "centroid". Detail see \code{\link[stats]{hclust}}.
#' @param font_size A numeric value specifying the font size for text labels in the dendrogram plots (default: \code{0.2}).
#' @param auto_cutree Logical; if \code{TRUE}, the function automatically determines the optimal number of clusters based on
#'     silhouette width (default: \code{FALSE}).
#' @param dend_title_size A numeric value specifying the font size of the dendrogram plot title (default: \code{15}).
#' @param cut_height_dist A numeric value used to adjust the vertical distance of the cut-off line annotation in the dendrogram
#'     plot (default: \code{0.2}).
#'
#' @return A list with two elements: \code{cluster_results} contains the hierarchical clustering objects for each group,
#'     and \code{cluster_figures} contains the corresponding ggplot2 dendrogram figures annotated with cluster information.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom cluster silhouette
#' @importFrom ggdendro dendro_data
#' @importFrom dendextend color_branches
#' @author Shijia Li
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming predicted_data is a list of predicted OTU abundance data frames:
#' cluster_output <- Data.cluster(predicted_data, clust_method = 'complete',
#'                                font_size = 0.3, auto_cutree = TRUE,
#'                                dend_title_size = 15, cut_height_dist = 0.2)
#' }
Data.cluster = function(predicted_data, clust_method, font_size = 0.2, auto_cutree = FALSE, dend_title_size = 15,cut_height_dist = 0.2) {

  cluster_results = list()
  cluster_figures = list()

  for (i in names(predicted_data)) {

    pred_data = as.data.frame(predicted_data[[i]][,-ncol(predicted_data[[i]])])

    if (is.null(pred_data) || !is.data.frame(pred_data) && !is.matrix(pred_data) || ncol(pred_data) < 3) {
      cluster_results[[i]] = NULL
      cat("Group",i,"was excluded due to having fewer than 3 OTUs, making it unsuitable for clustering.\n")
    } else {
      # 1. Calculate distance
      cor_matrix = cor(pred_data)
      pred_data_dist = as.dist((1-cor_matrix))

      # 2. Use hierachical cluster to do clustering of OTUs those have similar temporal curves
      pred_data_dist_hc = hclust(pred_data_dist, clust_method)

      if (auto_cutree == TRUE) {
        # 3. Cut dendrogram and calculate silhouette
        max_k = ncol(pred_data)-1
        sil_width = numeric(max_k)

        for (k in 2:max_k) {
          # Cut dendrogram
          cutree_k = cutree(pred_data_dist_hc,k)

          # Calculate silhouette
          sil = cluster::silhouette(cutree_k, pred_data_dist)
          sil_width[k] = mean(sil[,3])
        }

        # 4. Determine the optimal numbers of clusters
        optimal_k = which.max(sil_width)
        cat("Number of optimal clusters for group",i, ":", optimal_k,"\n")

        # 5. Visualise hcluster dendrogram with labeling of optimal clusters
        dend = as.dendrogram(pred_data_dist_hc)
        clusters = cutree(pred_data_dist_hc, optimal_k)[order.dendrogram(dend)]
        clusters.df = data.frame(label = names(clusters),cluster = factor(clusters))
        dend_color = dendextend::color_branches(dend = dend,clusters = clusters)
        dend_data <- ggdendro::dendro_data(dend,type = "rectangle")
        dend_data[["labels"]] <- merge(dend_data[["labels"]],clusters.df, by="label")

        cluster_figures[[i]] = ggplot()+
          geom_segment(data = dend_data$segments,aes(x = x,y = y,xend = xend,yend = yend))+
          geom_text(data = dend_data$labels,aes(x = x,y = y,label = label,color = cluster),
                    size = font_size,hjust = 0)+
          coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
          labs(title = paste("OTU profiles clustering plot for group ",i,"\n",
                             "Suggested optimal cluster: ", optimal_k, sep = ""))+
          theme(axis.line.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.x = element_blank(),
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title = element_text(size = dend_title_size),
                plot.margin = unit(c(0.1,1,0.1,0.1),'cm'))
      } else {
        dend = as.dendrogram(pred_data_dist_hc)
        dend_data <- ggdendro::dendro_data(dend,type = "rectangle")
        cluster_figures[[i]] = ggplot()+
          geom_segment(data = dend_data$segments,aes(x = x,y = y,xend = xend,yend = yend))+
          geom_text(data = dend_data$labels,aes(x = x,y = y,label = label),
                    size = font_size,hjust = 0)+
          coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
          labs(title = paste("OTU profiles clustering plot for group ",i,sep = ""))+
          theme(axis.line.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.x = element_blank(),
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title = element_text(size = dend_title_size),
                plot.margin = unit(c(0.1,1,0.1,0.1),'cm'))
        print(cluster_figures[[i]])
        cut_height = readline(paste("According to the dendrogram, the input cut-off line of ",i," is: ",sep = ""))
        cut_height = as.numeric(cut_height)
        cluster_figures_show = cluster_figures[[i]]+
          geom_hline(yintercept = cut_height,color = 'blue',linetype = 'dashed')+
          annotate("text",x = 3,y = cut_height+cut_height_dist,
                   label = paste('Cut-off line:',cut_height,sep = ""),color = 'blue')
        print(cluster_figures_show)
        cut_determine = readline(paste('Do you confirm the input cut-off value ',cut_height,' ? (Yes/No)',sep = ""))
        while (!cut_determine %in% c('Yes','YES','yes')) {
          cut_height = readline(paste("According to the dendrogram, the input cut-off line of ",i," is: ",sep = ""))
          cut_height = as.numeric(cut_height)
          cluster_figures_show = cluster_figures[[i]]+
            geom_hline(yintercept = cut_height,color = 'blue',linetype = 'dashed')+
            annotate("text",x = 3,y = cut_height+cut_height_dist,
                     label = paste('Cut-off line:',cut_height,sep = ""),color = 'blue')
          print(cluster_figures_show)
          cut_determine = readline(paste('Do you confirm the input cut-off value ',cut_height,' ? (Yes/No)',sep = ""))
        }
        cluster_figures[[i]] = cluster_figures_show
      }
      # 6. generate a list to store cluster result and visualization
      cluster_results[[i]] = pred_data_dist_hc
    }
  }
  # 7. return the clustered result list.
  cluster_data = list(cluster_results,cluster_figures)
  names(cluster_data) = c("cluster_results","cluster_figures")
  return(cluster_data)
}
