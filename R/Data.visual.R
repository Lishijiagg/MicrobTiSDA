#' @title Visualize Temporal OTU Profiles from Clustered Predicted Data
#' @description
#' The \code{Data.visual} function generates visualizations of temporal profiles for OTUs by integrating clustering results,
#' predicted time-series data, and design information. It produces ggplot2 figures for each group and for each cluster branch,
#' displaying smoothed curves of predicted OTU abundances over time. Optionally, the function overlays raw data points and fits linear
#' models to assess temporal trends, annotating the plots with model statistics when certain criteria are met.
#'
#' @details
#' This function uses hierarchical clustering results (obtained from a dendrogram) to cut the tree either by a specified height or
#'     by a user specified number of branches of each dendrogram in \code{cluster_results}. For each group in \code{cluster_results},
#'     the function extracts the corresponding predicted OTU data and raw design data. Temporal profiles are visualized by plotting
#'     smooth curves (using \code{stat_smooth}) for each cluster branch. When \code{plot_dots} is set to \code{TRUE}, the function
#'     overlays raw data points. Additionally, if \code{plot_lm} is \code{TRUE}, a linear model is fitted to the predicted data,
#'     and if the model meets specified thresholds for R-squared (\code{lm_R2}) and absolute slope (\code{lm_abs_slope})
#'     (i.e., R2 > 0.1 and absolute slope > 0.05), a dashed regression line is added along with an annotation of
#'     the R-squared and slope values. The resulting list of ggplot2 objects can be used to visually inspect the temporal dynamics of
#'     OTUs across different clusters and groups.
#'
#'
#' @param cluster_results A list object output from the \code{\link[MicrobTiSDA]{Data.cluster}}).
#' @param cutree_by A character string specifying the method to cut the dendrogram, either by \code{"height"} or by \code{"branches"}.
#' @param cluster_height A numeric vector specifying the cut-off height for each group when \code{cutree_by = "height"}.
#' @param cluster_branches A numeric vector specifying the number of clusters for each group when \code{cutree_by = "branches"}.
#' @param predicted_data The output data frame from the \code{\link[MicrobTiSDA]{Pred.data}}).
#' @param Design_data The output data from the \code{\link[MicrobTiSDA]{Design}}).
#' @param pre_processed_data The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
#'     pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.
#' @param Taxa A data frame providing taxonomic annotations for microbial species.
#' @param plot_dots Logical; if \code{TRUE}, raw data points are overlaid on the temporal curves (default: \code{TRUE}).
#' @param figure_x_scale A numeric value specifying the interval for x-axis breaks in the figures (default: \code{5}).
#' @param plot_lm Logical; if \code{TRUE}, a linear model is fitted to the predicted data to detect trends, and the regression line
#'     is added (default: \code{FALSE}).
#' @param lm_R2 A numeric threshold for the minimum R-squared value required to annotate the linear model (default: \code{0.01}).
#' @param lm_abs_slope A numeric threshold for the minimum absolute slope required to annotate the linear model (default: \code{0.005}).
#' @param title_size A numeric value specifying the font size for the plot title (default: \code{10}).
#' @param axis_title_size A numeric value specifying the font size for the axis titles (default: \code{8}).
#' @param axis_y_size A numeric value specifying the font size for the y-axis text (default: \code{5}).
#' @param axis_x_size A numeric value specifying the font size for the x-axis text (default: \code{5}).
#' @param lm_sig_size A numeric value specifying the font size for linear model annotation text (default: \code{5}).
#' @param legend_title_size A numeric value specifying the font size for legend titles (default: \code{5}).
#' @param legend_text_size A numeric value specifying the font size for legend text (default: \code{5}).
#' @param dots_size A numeric value specifying the size of the overlaid raw data points (default: \code{0.7}).
#' @importFrom ggplot2 ggplot
#' @return A list of lists of ggplot2 objects, where each top-level element corresponds to a subject and each sub-element corresponds
#'         to a cluster branch's temporal profile plot.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming cluster_results, predicted_data, Design_data, pre_processed_data, and Taxa are properly defined:
#' curves <- Data.visual(cluster_results = cluster_output,
#'                       cutree_by = "height",
#'                       cluster_height = c(0.3, 0.25, 0.4),
#'                       cluster_branches = NA,
#'                       predicted_data = predicted_data,
#'                       Design_data = design_data,
#'                       pre_processed_data = pre_processed_data,
#'                       Taxa = taxa_data,
#'                       plot_dots = TRUE,
#'                       figure_x_scale = 5,
#'                       plot_lm = TRUE,
#'                       lm_R2 = 0.01,
#'                       lm_abs_slope = 0.005,
#'                       title_size = 10,
#'                       axis_title_size = 8,
#'                       axis_y_size = 5,
#'                       axis_x_size = 5,
#'                       lm_sig_size = 5,
#'                       legend_title_size = 5,
#'                       legend_text_size = 5,
#'                       dots_size = 0.7)
#' # To view the plot for group 'Subject1' and cluster branch 2:
#' print(curves[['Subject1']][[2]])
#' }
Data.visual = function(cluster_results,
                      cutree_by = "height",
                      cluster_height = NA,
                      cluster_branches = NA,
                      predicted_data,
                      Design_data,
                      pre_processed_data,
                      Taxa = NULL,
                      plot_dots = TRUE,
                      figure_x_scale = 5,
                      plot_lm = FALSE,
                      lm_R2 = 0.01,
                      lm_abs_slope = 0.005,
                      title_size = 10,
                      axis_title_size = 8,
                      axis_y_size = 5,
                      axis_x_size = 5,
                      lm_sig_size = 5,
                      legend_title_size = 5,
                      legend_text_size = 5,
                      dots_size = 0.7) {

  curves_plot = list()
  cluster_results = cluster_results$cluster_results

  n = 0
  for (i in names(cluster_results)) {
    n = n + 1
    # 1. Select branches based on the dendrogram of clustering tree
    if (cutree_by == 'height') {
      cutrees = cutree(cluster_results[[i]], h = cluster_height[n])
    } else if (cutree_by == 'branches') {
      cutrees = cutree(cluster_results[[i]], k = cluster_branches[n])
    } else {
      stop("cutree_by must be either 'height' or 'branches'")
    }

    # 2. Extract the data from pre-processed data for plotting the transformed OTU count at each sampling time point
    rawdata = Design_data[Design_data[,i] == 1, ]

    ylim_min = min(min(pre_processed_data),
                   min(predicted_data[[i]][,-ncol(predicted_data[[i]])]))
    ylim_max = max(max(pre_processed_data),
                   max(predicted_data[[i]][,-ncol(predicted_data[[i]])]))

    # 3. Temporal profile visualizations
    curves_plot[[i]] = list()

    Pred_data_for_plot = predicted_data[[i]]
    for (j in as.numeric(levels(as.factor(cutrees)))) {
      cols_to_remove = names(which(cutrees != j))
      data_for_plot = Pred_data_for_plot[,!colnames(Pred_data_for_plot) %in% cols_to_remove]

      if(!is.null(Taxa)) {
        Annotate = Taxa[match(colnames(data_for_plot)[-ncol(data_for_plot)], rownames(Taxa)), ]
        # added mutate(across(where(is.factor), as.character)) %>%
        Annotate = Annotate %>% mutate(across(where(is.factor), as.character)) %>% mutate(across(everything(), ~na_if(.,"")))
        Taxon = apply(Annotate, 1, function(x) { tail(na.omit(x), 1) })
        Taxon = data.frame(Species_annotation = Taxon)
        colnames(data_for_plot) = c(paste(colnames(data_for_plot)[-ncol(data_for_plot)],
                                          " (", Taxon[, 1], ")", sep = ""),
                                    "Predicted_Time")
      } else {
        colnames(data_for_plot) = colnames(data_for_plot)
      }

      data_for_plot_long = data_for_plot %>%
        pivot_longer(cols = -Predicted_Time,names_to = "Features",values_to = "Abundance")
      curves_plot[[i]][[j]] = ggplot(data_for_plot_long, aes(x = Predicted_Time, y = Abundance, color = Features)) +
        stat_smooth(se = FALSE) +
        labs(title = paste(i,"_cluster_",j,sep = ""),
             x = "Time",y = "MCLR transformed Abundance")+
        theme_bw()+
        theme(title = element_text(size = title_size),
              axis.title = element_text(size = axis_title_size),
              axis.text.y = element_text(size = axis_y_size),
              axis.text.x = element_text(size = axis_x_size),
              legend.title = element_text(size = legend_title_size),
              legend.text = element_text(size = legend_text_size))+
        scale_x_continuous(breaks = seq(1,max(data_for_plot_long$Predicted_Time),figure_x_scale),
                           labels = seq(1,max(data_for_plot_long$Predicted_Time),figure_x_scale))+
        ylim(ylim_min-2,ylim_max+2)

      if (plot_dots == TRUE) {
        cols_to_reserve = names(which(cutrees == j))
        point_data_for_plot = rawdata[colnames(rawdata) %in% cols_to_reserve]
        point_data_for_plot$Time = rawdata$Time

        if (!is.null(Taxa)) {
          Annotate = Taxa[match(colnames(point_data_for_plot)[-ncol(point_data_for_plot)],
                                rownames(Taxa)),]
          Annotate = Annotate %>% mutate(across(where(is.factor), as.character)) %>% mutate(across(everything(), ~na_if(.,"")))
          Taxon = apply(Annotate,1,function(x) {tail(na.omit(x),1)})
          Taxon = data.frame(Species_annotation = Taxon)
          colnames(point_data_for_plot) = c(paste(colnames(point_data_for_plot)[-ncol(point_data_for_plot)],
                                                  " (",Taxon[,1],")",sep = ""),
                                            "Time")
        } else {
          colnames(point_data_for_plot) = colnames(point_data_for_plot)
        }

        point_data_for_plot_long = point_data_for_plot %>%
          pivot_longer(cols = -Time,names_to = "Features",values_to = "Abundance")
        curves_plot[[i]][[j]] = curves_plot[[i]][[j]]+
          geom_point(data = point_data_for_plot_long,
                     aes(x = Time,
                         y = Abundance,
                         color = Features),size = dots_size)

      } else {
        curves_plot[[i]][[j]] = curves_plot[[i]][[j]]
      }

      if (plot_lm == TRUE) {
        lm_model = lm(Abundance ~ Predicted_Time,data = data_for_plot_long)
        lm_model_summary = summary(lm_model)
        r_squared = lm_model_summary$r.squared
        slope = as.numeric(lm_model$coefficients[2])
        p_value = lm_model_summary$coefficients[, "Pr(>|t|)"]
        slope_p = as.numeric(p_value[2])
        if (slope_p < 0.05 && r_squared > lm_R2 && abs(slope) > lm_abs_slope) {
          if (slope > 0) {
            curves_plot[[i]][[j]] = curves_plot[[i]][[j]]+
              stat_smooth(data = data_for_plot_long,aes(x = Predicted_Time,y = Abundance),
                          method = 'lm',linetype = 'dashed',se = TRUE,color = 'darkgrey')+
              annotate("text",x=-Inf,y=Inf,label = paste("Slope=",round(slope,4)),
                       hjust = -1.1, vjust = 1.5,size = lm_sig_size,color = 'darkgray')
          } else {
            curves_plot[[i]][[j]] = curves_plot[[i]][[j]]+
              stat_smooth(data = data_for_plot_long,aes(x = Predicted_Time,y = Abundance),
                          method = 'lm',linetype = 'dashed',se = TRUE,color = 'darkgrey')+
              annotate("text",x=-Inf,y=Inf,label = paste("Slope=",round(slope,4)),
                       hjust = -1.1, vjust = 1.5,size = lm_sig_size,color = 'darkgray')
          }
        } else {
          curves_plot[[i]][[j]] = NULL
        }
      }
    }
  }

  # 4. Return the visualized temporal profiles of OTUs for each group
  return(curves_plot)
}
