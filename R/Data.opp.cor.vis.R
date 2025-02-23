#' @title Visualize Opposing Temporal profiles among Microbial features.
#' @description
#' This function identifies and visualizes OTUs/ASVs that exhibit opposing temporal trends based on a correlation threshold.
#'     It computes the correlation matrix of predicted OTU/ASV time-series data, selects those OTUs with correlations below a specified threshold,
#'     and generates smoothed temporal profile plots. Optionally, raw data points are overlaid and taxonomic annotations are added if provided.
#'
#' @details
#' For each group in the \code{predicted_data} list, the function first removes the time column and computes a correlation matrix
#' from the predicted OTU data. It then extracts, for each OTU, the subset of OTUs that show a correlation lower than the specified threshold
#' (\eqn{\text{ng\_cor\_thres}}). If an OTU does not have any opposing trends (i.e., all correlations exceed the threshold), it is skipped.
#' For those OTUs meeting the criteria, the function restructures the data into long format and plots the temporal profiles using
#' \code{geom_smooth} to display the primary trend (solid line) and opposing trends (dashed lines). If \code{plot_dots} is \code{TRUE},
#' the raw data points extracted from the design data are also overlaid. When taxonomic annotations are provided via \code{Taxa},
#' OTU labels are augmented with species information. The x-axis is scaled according to \code{figure_x_scale}, and plot aesthetics
#' (titles, axis text, legends, and dot sizes) can be customized using the respective parameters.
#'
#' @param predicted_data The output data frame from the \code{\link[MicrobTiSDA]{Pred.data}}).
#' @param pre_processed_data The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
#'     pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.
#' @param Design_data The output data from the \code{\link[MicrobTiSDA]{Design}}).
#' @param ng_cor_thres A numeric value specifying the correlation threshold below which OTUs are considered to exhibit opposing trends.
#' @param Taxa A data frame providing taxonomic annotations for microbial species.
#' @param plot_dots Logical; if \code{TRUE}, raw data points are overlaid on the temporal curves (default: \code{TRUE}).
#' @param figure_x_scale A numeric value specifying the interval for x-axis breaks in the figures (default: \code{5}).
#' @param title_size A numeric value specifying the font size for the plot title (default: \code{10}).
#' @param axis_title_size A numeric value specifying the font size for the axis titles (default: \code{8}).
#' @param axis_text_y_size A numeric value specifying the font size for the y-axis text (default: \code{5}).
#' @param axis_text_x_size A numeric value specifying the font size for the x-axis text (default: \code{5}).
#' @param legend_title_size A numeric value specifying the font size for legend titles (default: \code{5}).
#' @param legend_text_size A numeric value specifying the font size for legend text (default: \code{5}).
#' @param dots_size A numeric value specifying the size of the overlaid raw data points (default: \code{0.7}).
#'
#' @author Shijia Li
#' @return A list of lists of ggplot2 objects. The top-level list is keyed by group names, and each sublist contains plots for each OTU
#'         that exhibits opposing trends (i.e., has correlations below \code{ng_cor_thres}) in that group.
#' @export
#'
Data.opp.cor.vis = function(predicted_data,
                            pre_processed_data,
                            Design_data,
                            ng_cor_thres,
                            Taxa = NULL,
                            plot_dots = TRUE,
                            figure_x_scale = 5,
                            title_size = 10,
                            axis_title_size = 10,
                            axis_text_y_size = 8,
                            axis_text_x_size = 8,
                            legend_title_size = 10,
                            legend_text_size = 8,
                            dots_size = 0.6){

  curves_plot = list()

  for (i in names(predicted_data)) {
    pred_data = as.data.frame(predicted_data[[i]])
    rawdata = Design_data[Design_data[,i] == 1, ]
    curves_plot[[i]] = list()

    ylim_min = min(min(pre_processed_data),
                   min(predicted_data[[i]][,-ncol(predicted_data[[i]])]))
    ylim_max = max(max(pre_processed_data),
                   max(predicted_data[[i]][,-ncol(predicted_data[[i]])]))

    pred_data = pred_data[,-ncol(pred_data)]
    cor_matrix = cor(pred_data)
    oppo_cor = as.data.frame(cor_matrix)

    for (j in colnames(oppo_cor)) {
      oppo_cor_otu = oppo_cor[j]
      if (min(oppo_cor_otu) > ng_cor_thres) {
        curves_plot[[i]][[j]] = NULL
        next
      } else {
        oppo_otu = subset(oppo_cor_otu,oppo_cor_otu[j] < ng_cor_thres)
        oppo_otu_vec = c(colnames(oppo_otu),rownames(oppo_otu))

        data_for_plot = pred_data[,oppo_otu_vec]
        data_for_plot$Predicted_Time = predicted_data[[i]]$Predicted_Time

        if(!is.null(Taxa)) {
          Annotate = Taxa[match(colnames(data_for_plot)[-ncol(data_for_plot)], rownames(Taxa)), ]
          # added mutate(across(where(is.factor), as.character)) %>%
          Annotate = Annotate %>% mutate(across(where(is.factor), as.character)) %>% mutate(across(everything(), ~na_if(.,"")))
          Taxon = apply(Annotate, 1, function(x) { tail(na.omit(x), 1) })
          Taxon = data.frame(Species_annotation = Taxon)
          colnames(data_for_plot) = c(paste(colnames(data_for_plot)[-ncol(data_for_plot)],
                                            " (", Taxon[, 1], ")", sep = ""),
                                      "Predicted_Time")
        }

        data_for_plot_long = data_for_plot %>%
          pivot_longer(cols = -Predicted_Time,names_to = "Features",values_to = "Abundance")

        data_for_plot_long_1 = subset(data_for_plot_long,data_for_plot_long$Features == colnames(data_for_plot[1]))
        data_for_plog_long_2 = subset(data_for_plot_long,data_for_plot_long$Features != colnames(data_for_plot[1]))

        curves_plot[[i]][[j]] = ggplot(data_for_plot_long_1, aes(x = Predicted_Time, y = Abundance, color = Features)) +
          geom_smooth(se = FALSE) +
          geom_smooth(data = data_for_plog_long_2, aes(x = Predicted_Time, y = Abundance, color = Features),
                      se = FALSE, linetype = 'dashed') +
          labs(title = paste('Features in ', i, ' with opposite trends to ', j, sep = ""),
               x = "Time", y = "Predicted values") +
          theme_bw()+
          theme(plot.title = element_text(size = title_size),
                axis.title = element_text(size = axis_title_size),
                axis.text.y = element_text(size = axis_text_y_size),
                axis.text.x = element_text(size = axis_text_x_size),
                legend.title = element_text(size = legend_title_size),
                legend.text = element_text(size = legend_text_size)) +
          scale_x_continuous(breaks = seq(1, max(data_for_plot_long$Predicted_Time), figure_x_scale),
                             labels = seq(1, max(data_for_plot_long$Predicted_Time), figure_x_scale)) +
          ylim(ylim_min - 2, ylim_max + 2)

        if (plot_dots == TRUE) {
          point_data_for_plot = rawdata[colnames(rawdata) %in% oppo_otu_vec]
          point_data_for_plot$Time = rawdata$Time

          if(!is.null(Taxa)) {
            Annotate = Taxa[match(colnames(point_data_for_plot)[-ncol(point_data_for_plot)], rownames(Taxa)), ]
            # added mutate(across(where(is.factor), as.character)) %>%
            Annotate = Annotate %>% mutate(across(where(is.factor), as.character)) %>% mutate(across(everything(), ~na_if(.,"")))
            Taxon = apply(Annotate, 1, function(x) { tail(na.omit(x), 1) })
            Taxon = data.frame(Species_annotation = Taxon)
            colnames(point_data_for_plot) = c(paste(colnames(point_data_for_plot)[-ncol(point_data_for_plot)],
                                                    " (", Taxon[, 1], ")", sep = ""),
                                              "Time")
          }
          point_data_for_plot_long = point_data_for_plot %>%
            pivot_longer(cols = -Time,names_to = "Features",values_to = "Abundance")
          curves_plot[[i]][[j]] = curves_plot[[i]][[j]]+
            geom_point(data = point_data_for_plot_long,
                       aes(x = Time,
                           y = Abundance,
                           color = Features),size = dots_size)
        }
      }
    }

  }
  return(curves_plot)
}
