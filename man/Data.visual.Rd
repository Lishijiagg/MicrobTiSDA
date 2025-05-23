% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.visual.R
\name{Data.visual}
\alias{Data.visual}
\title{Visualize Temporal OTU Profiles from Clustered Predicted Data}
\usage{
Data.visual(
  cluster_results,
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
  dots_size = 0.7
)
}
\arguments{
\item{cluster_results}{A list object output from the \code{\link[MicrobTiSDA]{Data.cluster}}).}

\item{cutree_by}{A character string specifying the method to cut the dendrogram, either by \code{"height"} or by \code{"branches"}.}

\item{cluster_height}{A numeric vector specifying the cut-off height for each group when \code{cutree_by = "height"}.}

\item{cluster_branches}{A numeric vector specifying the number of clusters for each group when \code{cutree_by = "branches"}.}

\item{predicted_data}{The output data frame from the \code{\link[MicrobTiSDA]{Pred.data}}).}

\item{Design_data}{The output data from the \code{\link[MicrobTiSDA]{Design}}).}

\item{pre_processed_data}{The transformed data output from the \code{\link[MicrobTiSDA]{Data.trans}} function. A
pre-processed OTU data frame with sample IDs as row names and OTU IDs as column names.}

\item{Taxa}{A data frame providing taxonomic annotations for microbial species.}

\item{plot_dots}{Logical; if \code{TRUE}, raw data points are overlaid on the temporal curves (default: \code{TRUE}).}

\item{figure_x_scale}{A numeric value specifying the interval for x-axis breaks in the figures (default: \code{5}).}

\item{plot_lm}{Logical; if \code{TRUE}, a linear model is fitted to the predicted data to detect trends, and the regression line
is added (default: \code{FALSE}).}

\item{lm_R2}{A numeric threshold for the minimum R-squared value required to annotate the linear model (default: \code{0.01}).}

\item{lm_abs_slope}{A numeric threshold for the minimum absolute slope required to annotate the linear model (default: \code{0.005}).}

\item{title_size}{A numeric value specifying the font size for the plot title (default: \code{10}).}

\item{axis_title_size}{A numeric value specifying the font size for the axis titles (default: \code{8}).}

\item{axis_y_size}{A numeric value specifying the font size for the y-axis text (default: \code{5}).}

\item{axis_x_size}{A numeric value specifying the font size for the x-axis text (default: \code{5}).}

\item{lm_sig_size}{A numeric value specifying the font size for linear model annotation text (default: \code{5}).}

\item{legend_title_size}{A numeric value specifying the font size for legend titles (default: \code{5}).}

\item{legend_text_size}{A numeric value specifying the font size for legend text (default: \code{5}).}

\item{dots_size}{A numeric value specifying the size of the overlaid raw data points (default: \code{0.7}).}
}
\value{
A list of lists of ggplot2 objects, where each top-level element corresponds to a subject and each sub-element corresponds
to a cluster branch's temporal profile plot.
}
\description{
The \code{Data.visual} function generates visualizations of temporal profiles for OTUs by integrating clustering results,
predicted time-series data, and design information. It produces ggplot2 figures for each group and for each cluster branch,
displaying smoothed curves of predicted OTU abundances over time. Optionally, the function overlays raw data points and fits linear
models to assess temporal trends, annotating the plots with model statistics when certain criteria are met.
}
\details{
This function uses hierarchical clustering results (obtained from a dendrogram) to cut the tree either by a specified height or
by a user specified number of branches of each dendrogram in \code{cluster_results}. For each group in \code{cluster_results},
the function extracts the corresponding predicted OTU data and raw design data. Temporal profiles are visualized by plotting
smooth curves (using \code{stat_smooth}) for each cluster branch. When \code{plot_dots} is set to \code{TRUE}, the function
overlays raw data points. Additionally, if \code{plot_lm} is \code{TRUE}, a linear model is fitted to the predicted data,
and if the model meets specified thresholds for R-squared (\code{lm_R2}) and absolute slope (\code{lm_abs_slope})
(i.e., R2 > 0.1 and absolute slope > 0.05), a dashed regression line is added along with an annotation of
the R-squared and slope values. The resulting list of ggplot2 objects can be used to visually inspect the temporal dynamics of
OTUs across different clusters and groups.
}
\examples{
\dontrun{
# Assuming cluster_results, predicted_data, Design_data, pre_processed_data, and Taxa are properly defined:
curves <- Data.visual(cluster_results = cluster_output,
                      cutree_by = "height",
                      cluster_height = c(0.3, 0.25, 0.4),
                      cluster_branches = NA,
                      predicted_data = predicted_data,
                      Design_data = design_data,
                      pre_processed_data = pre_processed_data,
                      Taxa = taxa_data,
                      plot_dots = TRUE,
                      figure_x_scale = 5,
                      plot_lm = TRUE,
                      lm_R2 = 0.01,
                      lm_abs_slope = 0.005,
                      title_size = 10,
                      axis_title_size = 8,
                      axis_y_size = 5,
                      axis_x_size = 5,
                      lm_sig_size = 5,
                      legend_title_size = 5,
                      legend_text_size = 5,
                      dots_size = 0.7)
# To view the plot for group 'Subject1' and cluster branch 2:
print(curves[['Subject1']][[2]])
}
}
