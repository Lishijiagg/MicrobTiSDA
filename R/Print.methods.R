#' Print Method for TransformedData Object
#'
#' @param x A \code{TransformedData} object
#' @param ... Additional arguments
#'
#' @method print TransformedData
#' @export
print.TransformedData <- function(x, ...) {
  cat("An object of class 'TransformedData'\n")
  cat("Data dimensions: ", dim(x$transformed_data)[1], " x ", dim(x$transformed_data)[2], "\n", sep = "")
}


#' Print Method for FilteredData Object
#'
#' @param x A \code{FilteredData} object
#' @param ... Additional arguments
#'
#' @method print FilteredData
#' @export
print.FilteredData <- function(x, ...) {
  cat("FilteredData object\n")
  cat("---------------------------------\n")
  cat("Number of OTUs: ", nrow(x$filtered_table), "\n")
  cat("Number of Samples: ", ncol(x$filtered_table), "\n")
  cat("Parameters used:\n")
  print(x$parameters)
  invisible(x)
}


#' Print MicrobTiSDA.interpolate object
#'
#' @param x An object of class \code{MicrobTiSDA.interpolate}
#' @param ... Additional arguments
#' @method print MicrobTiSDA.interpolate
#'
#' @export
print.MicrobTiSDA.interpolate <- function(x, ...) {
  cat("Interpolated microbiome time-series data\n")
  cat("Interpolation method:", x$Method, "\n")
  cat("Number of features:", nrow(x$Interpolated_Data), "\n")
  cat("Number of interpolated samples:", ncol(x$Interpolated_Data), "\n")
  invisible(x)
}


#' Print MicrobTiSDA objects
#'
#' @param x An object of class \code{microbTiSDA}
#' @param ... Additional arguments
#' @method print microbTiSDA
#' @export
print.microbTiSDA <- function(x, ...) {
  cat("MicrobTiSDA Object of class:", class(x), "\n")
  if (!is.null(x$params)) {
    cat("Parameters used:\n")
    print(x$params)
  }
  if (!is.null(x$results)) {
    cat("Number of groups/items in results:", length(x$results), "\n")
    cat("Use summary() for a detailed overview.\n")
  }
  invisible(x)
}


#' Print information of species interaction dynamic visualizations
#'
#' @param x An object of class microbTiSDA_dynamic_vis
#' @param ... Additional arguments
#'
#' @method print microbTiSDA_dynamic_vis
#' @export
print.microbTiSDA_dynamic_vis <- function(x, ...) {
  cat("MicrobTiSDA Dynamic Interaction Visualization\n")
  cat("Number of groups: ", length(x$visualization), "\n")
  cat("Threshold: ", x$params$threshold, "\n")
  cat("Core arrow number: ", x$params$core_arrow_num, "\n")
  cat("Font size: ",x$params$fontsize,"\n")
  invisible(x)
}


#' Print information of the Design matrix
#'
#' @param x An object of class Design.
#' @param ... Additional arguments
#' @method print Design
#' @export
print.Design <- function(x, ...) {
  cat("Design object\n")
  cat("Number of samples:", nrow(x$data), "\n")
  cat("Number of predictors:", ncol(x$data) - 2, " (excluding Time and ID)\n") # exclude Time & ID
  if (!is.null(x$params$Group_var)) {
    cat("Grouping variable(s):", paste(x$params$Group_var, collapse = ", "), "\n")
  } else {
    cat("No grouping variable specified\n")
  }
  invisible(x)
}


#' Print information of fitted natural spline regression models
#'
#' @param x An object of class MicrobTiSDA_spline_regression
#' @param ... Additional arguments
#' @method print MicrobTiSDA_spline_regression
#' @export
print.MicrobTiSDA_spline_regression <- function(x, ...) {
  cat("MicrobTiSDA Spline Regression Object\n")
  cat("Number of independent variables:", length(x$fitted_model), "\n")
  cat("Knots info available for each model\n")
  invisible(x)
}


#' Print the information of fitted regression model predicted data
#'
#' @param x An object of class PredictedData
#' @param ... Additional arguments
#'
#' @method print PredictedData
#' @export
print.PredictedData <- function(x, ...) {
  cat("PredictedData object\n")
  cat("Number of groups:", length(x), "\n")
  cat("Groups:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
}


#' Print the information of fitted temporal profile clusters
#'
#' @param x An object of class MicrobTiSDA.cluster
#' @param ... Additional arguments
#' @method print MicrobTiSDA.cluster
#' @export
print.MicrobTiSDA.cluster <- function(x, ...) {
  cat("MicrobTiSDA Clustering Object\n")
  cat("Groups clustered:", length(x$cluster_results), "\n")
  cat("Use summary() to see details, or plot() to visualize.\n")
  invisible(x)
}


#' Print the plot information of user selected feature clustering figures
#'
#' @param x An object of class MicrobTiSDA.clusterCut
#' @param ... Additional arguments
#' @method print MicrobTiSDA.clusterCut
#' @export
print.MicrobTiSDA.clusterCut <- function(x, ...) {
  cat("MicrobTiSDA Cluster Cutting Object\n")
  cat("Groups:", length(x$cluster_results), "\n")
  cat("Use summary() to check number of clusters per group, or plot() to visualize.\n")
  invisible(x)
}

#' Print the information of those clustered microbial features' temporal patterns
#'
#' @param x An object of class MicrobTiSDA.visual
#' @param ... Additional arguments
#' @method print MicrobTiSDA.visual
#' @export
print.MicrobTiSDA.visual <- function(x, ...) {
  cat("MicrobTiSDA Visualization Object\n")
  cat("Groups:", length(x$curves_plot), "\n")
  cat("Use summary() to check available clusters, or plot() to render figures.\n")
  invisible(x)
}


#' Information of Microbial features with opposite temporal shapes to users selected feature
#'
#' @param x An object of class DataOppCorVis
#' @param ... Additional arguments
#' @method print DataOppCorVis
#' @export
print.DataOppCorVis <- function(x, ...) {
  cat("S3 Object of class 'DataOppCorVis'\n")
  cat("Groups:", length(x$curves_plot), "\n")
  invisible(x)
}

#' Information of the constructed Random Forest classification model
#'
#' @param x An object of class DataRFClassifier
#' @param ... Additional arguments
#' @method print DataRFClassifier
#' @export
print.DataRFClassifier <- function(x, ...) {
  message("S3 Object of class 'DataRFClassifier'")
  message("Contains the trained random forest model, feature importance, and cross-validation results.")
  invisible(x)
}

#' Information of the fitted mixed-effect spline regression models
#'
#' @param x An object of class RegMESR
#' @param ... Additional arguments
#' @method print RegMESR
#' @export
print.RegMESR <- function(x, ...) {
  cat("S3 Object of class 'RegMESR'\n")
  cat("Contains:\n")
  cat(" - Fitted GAM models for each OTU\n")
  cat(" - Knot information for each model\n")
  cat(" - Model parameters\n")
  invisible(x)
}


#' Print method for PredictedDataMESR objects
#'
#' @param x An object of class PredictedDataMESR.
#' @param ... Additional arguments.
#'
#' @method print PredictedDataMESR
#' @export
print.PredictedDataMESR <- function(x, ...) {
  cat("PredictedDataMESR object\n")
  cat("Number of groups:", length(x$groups), "\n")
  cat("Groups:", paste(x$groups, collapse = ", "), "\n")
  cat("Time variable:", x$params$Sample_Time, "\n")
  cat("Prediction stored for each group in `object$data`\n")
  invisible(x)
}


#' Information of visualizations of feature temporal patterns fitted with mixed-effect spline regression models
#'
#' @param x An object of class MicrobTiSDA.MESRvisual
#' @param ... Additional arguments
#' @method print MicrobTiSDA.MSERvisual
#' @export
print.MicrobTiSDA.MSERvisual <- function(x, ...) {
  cat("S3 Object of class 'MicrobTiSDA.MSERvisual'\n")
  cat("Contains visualization results of MESR clustering.\n")
  cat("Number of groups:", length(x$plots), "\n")
  invisible(x)
}
