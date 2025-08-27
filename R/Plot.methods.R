#' Plot Method for microbTiSDA Dynamic Visualization Objects
#'
#' This function iterates over all visualization elements stored in a
#' \code{microbTiSDA_dynamic_vis} object and prints each plot with its name.
#'
#' @param x A \code{microbTiSDA_dynamic_vis} object containing plots in the
#'   \code{visualization} list element.
#' @param groups A character vector specifying the names of groups to plot.
#'     if \code{NULL}, all groups in \code{x$visualization} will be plotted.
#' @param ... Additional arguments
#'
#' @export
plot.microbTiSDA_dynamic_vis <- function(x, groups = NULL, ...) {
  if (is.null(groups)) {
    groups <- names(x$visualization)
  }
  for (g in groups) {
    cat("Plot:", names(x$visualization)[g],"\n")
    print(x$visualization[[g]])
  }
  invisible(x)
}

#' Plot Method for MicrobTiSDA Cluster Objects
#'
#' This function prints cluster plots stored in a \code{MicrobTiSDA.cluster} object.
#' Users can specify which groups to plot; if \code{groups} is \code{NULL},
#' all available cluster plots will be displayed.
#' @param x A \code{MicrobTiSDA.cluster} object containing cluster plots
#'   in the \code{cluster_figures} list element.
#' @param groups A character vector specifying the names of groups to plot.
#'   If \code{NULL}, all groups in \code{x$cluster_figures} will be plotted.
#' @param ... Additional arguments (currently not used) for compatibility with
#'   generic \code{plot} methods.
#'
#' @export
plot.MicrobTiSDA.cluster <- function(x, groups = NULL, ...) {
  if (is.null(groups)) {
    groups <- names(x$cluster_figures)
  }
  for (g in groups) {
    if (!is.null(x$cluster_figures[[g]])) {
      print(x$cluster_figures[[g]])
    }
  }
  invisible(x)
}


#' Plot Method for MicrobTiSDA Cluster Cut Objects
#'
#' This function prints cluster plots stored in a \code{MicrobTiSDA.clusterCut} object.
#' Users can specify which groups to plot; if \code{groups} is \code{NULL},
#' all available cluster plots will be displayed.
#' @param x x A \code{MicrobTiSDA.clusterCut} object containing cluster plots
#'   in the \code{cluster_figures} list element.
#' @param groups A character vector specifying the names of groups to plot.
#'   If \code{NULL}, all groups in \code{x$cluster_figures} will be plotted.
#' @param ... Additional arguments (currently not used) for compatibility with
#'   generic \code{plot} methods.
#'
#' @export
plot.MicrobTiSDA.clusterCut <- function(x, groups = NULL, ...) {
  if (is.null(groups)) groups <- names(x$cluster_figures)
  for (g in groups) {
    if (!is.null(x$cluster_figures[[g]])) {
      print(x$cluster_figures[[g]])
    }
  }
  invisible(x)
}

#' Plot Method for MicrobTiSDA Visual Objects
#'
#' This function prints feature cluster plots stored in a \code{MicrobTiSDA.visual} object.
#' Users can specify which groups and which clusters within each group to plot.
#' If \code{groups} is \code{NULL}, all groups will be displayed. If \code{clusters}
#' is \code{NULL}, all clusters within the selected groups will be plotted.
#' @param x A \code{MicrobTiSDA.visual} object containing feature cluster plots
#'   in the \code{plots} list element.
#' @param groups A character vector specifying the names of groups to plot.
#'   If \code{NULL}, all groups in \code{x$plots} will be plotted.
#' @param clusters An integer vector specifying which clusters to plot within
#'   each selected group. If \code{NULL}, all clusters are plotted.
#' @param ... Additional arguments (currently not used) for compatibility with
#'   generic \code{plot} methods.
#'
#' @export
plot.MicrobTiSDA.visual <- function(x, groups = NULL, clusters = NULL, ...) {
  if (is.null(groups)) groups <- names(x$plots)
  for (g in groups) {
    cat(paste0("In ", g, ", number of feature clusters -- ",length(x$plots[[g]]),"\n"))
    if (!is.null(x$plots[[g]])) {
      these_clusters <- if (is.null(clusters)) c(1:length(x$plots[[g]])) else clusters
      for (c in these_clusters) {
        if (!is.null(x$plots[[g]][[c]])) {
          print(x$plots[[g]][[c]])
        }
      }
    }
  }
  invisible(x)
}

#' Plot Method for DataOppCorVis Objects
#'
#' This function prints correlation curve plots stored in a \code{DataOppCorVis} object.
#' Users can specify a group and optionally a specific feature within that group to plot.
#' @param x A \code{DataOppCorVis} object containing correlation curve plots
#'   in the \code{curves_plot} list element.
#' @param group A character string specifying the group to plot. Must be provided.
#' @param feature A character string specifying the feature within the group to plot.
#'   If omitted, all features within the specified group will be plotted.
#' @param ... Additional arguments
#'
#' @export
plot.DataOppCorVis <- function(x, group, feature, ...) {
  if (!missing(group) && !missing(feature)) {
    print(x$curves_plot[[group]][[feature]])
  } else if (!missing(group)) {
    for (fig in names(x$curves_plot[[group]])) {
      print(x$curves_plot[[group]][[fig]])
    }
  }
  else {
    stop("Specify group and feature")
  }
}

#' Plot Method for DataRFClassifier Objects
#'
#' This function plots the margin scores and cross-validation curve of a
#' \code{DataRFClassifier} object.
#'
#' @param x A \code{DataRFClassifier} object containing margin scores and
#'   cross-validation curve plots.
#' @param ... Additional arguments
#'
#' @export
plot.DataRFClassifier <- function(x,...) {
  if (!inherits(x, "DataRFClassifier")) stop("Input must be a DataRFClassifier object.")

  p1 <- NULL
  p2 <- NULL

  # Margin scores plot
  if (!is.null(x$Margin_scores_train)) {
    p1 <- x$Margin_scores_train + ggtitle("Margin Scores (Train Set)")
  } else {
    warning("Margin scores plot not found in the object.")
  }

  # Cross-validation curve plot
  if (!is.null(x$cross_validation)) {
    p2 <- x$cross_validation + ggtitle("Cross-Validation Curve")
  } else {
    warning("Cross-validation plot not found in the object.")
  }

  figures <- list(p1,p2)
  names(figures) <- c("Margin Scores","CV-curve")

  print(figures)

  invisible(figures)
}


#' Plot Method for RfBiomarker Objects
#'
#' This function prints the cross-validation figure of a \code{RfBiomarker} object
#' and displays basic information about the object contents.
#'
#' @param x A \code{RfBiomarker} object containing selected important OTUs,
#'   predictions on training and test sets, and a cross-validation figure.
#' @param ... Additional arguments
#'
#' @export
plot.RfBiomarker <- function(x, ...) {
  print(x$cross_validation_fig)
  invisible(x)
}

#' Plot Method for MicrobTiSDA MSER Visualization Objects
#'
#' This function prints feature cluster plots stored in a
#' \code{MicrobTiSDA.MSERvisual} object. Users can specify which groups and
#' which clusters within each group to plot. If \code{groups} or \code{clusters}
#' is \code{NULL}, all groups or clusters are plotted.
#' @param x A \code{MicrobTiSDA.MSERvisual} object containing feature cluster plots
#'   in the \code{plots} list element.
#' @param groups A character vector specifying the names of groups to plot.
#'   If \code{NULL}, all groups in \code{x$plots} will be plotted.
#' @param clusters An integer vector specifying which clusters to plot within
#'   each selected group. If \code{NULL}, all clusters are plotted.
#' @param ... Additional arguments
#'
#' @export
plot.MicrobTiSDA.MSERvisual <- function(x, groups = NULL, clusters = NULL, ...) {
  if (is.null(groups)) groups <- names(x$plots)
  for (g in groups) {
    if (!is.null(x$plots[[g]])) {
      these_clusters <- if (is.null(clusters)) c(1:length(x$plots[[g]])) else clusters
      for (c in these_clusters) {
        if (!is.null(x$plots[[g]][[c]])) {
          print(x$plots[[g]][[c]])
        }
      }
    }
  }
  invisible(x)
}
