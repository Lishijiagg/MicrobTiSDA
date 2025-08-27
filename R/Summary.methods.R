#' Summary Method for FilteredData Object
#'
#' @param object A \code{FilteredData} object.
#' @param ... Additional arguments
#'
#' @return A summary list containing dimensions and filtering parameters
#' @method summary FilteredData
#' @export
summary.FilteredData <- function(object, ...) {
  summary_list <- list(
    n_OTUs = nrow(object$filtered_table),
    n_Samples = ncol(object$filtered_table),
    parameters = object$parameters
  )
  cat("Summary of FilteredData object\n")
  cat("---------------------------------\n")
  cat("Number of OTUs: ", summary_list$n_OTUs, "\n")
  cat("Number of Samples: ", summary_list$n_Samples, "\n")
  cat("Filtering Parameters:\n")
  print(summary_list$parameters)
  cat("Preview: \n")
  print(head(object$filtered_table))
  invisible(summary_list)
}


#' Summary method for MicrobTiSDA.interpolate objects
#'
#' Provides a structure summary of a \code{MicrobTiSDA.interpolate} object, including its parameters and results.
#' It prints object class, the information of interpolated microbial count data and updated metadata.
#'
#' @param object A \code{MicrobTiSDA.interpolate} object created by the package functions.
#' @param ... Additional arguments
#'
#' @export
summary.MicrobTiSDA.interpolate <- function(object, ...) {
  cat("Summary of MicrobTiSDA.interpolate object\n")
  cat("------------------------------------------\n")
  cat("Interpolation method:", object$Method, "\n")
  cat("Features:", nrow(object$Interpolated_Data), "\n")
  cat("Samples:", ncol(object$Interpolated_Data), "\n")
  cat("Preview interpolated data: \n")
  print(head(object$Interpolated_Data))
  cat("------------------------------------------\n")
  cat("Metadata rows:", nrow(object$Interpolated_Data_metadata), "\n")
  cat("Preview interpolated metadata: \n")
  print(head(object$Interpolated_Data_metadata))
  cat("Groups:", length(unique(object$Interpolated_Data_metadata$Group)), "\n")
  invisible(object)
}


#' Summary Method for TransformedData Objects
#'
#' @param object A \code{TransformedData} object returned by \code{Data.trans()}.
#' @param ... Additional arguments.
#'
#' @export
summary.TransformedData <- function(object, ...) {
  cat("Summary of TransformedData object\n")
  cat("---------------------------------\n")
  cat("Number of samples (rows):", nrow(object$transformed_data), "\n")
  cat("Number of features (columns):", ncol(object$transformed_data), "\n")

  cat("Preview of transformed data:\n")
  print(object$transformed_data[1:10,1:10])


  invisible(object)
}


#' Summary method for MicrobTiSDA objects
#' Provides a structure summary of a \code{microbTiSDA} object, including its parameters and results. It prints object class,
#' parameters settings, number of result groups/items, and previews of tabular of list-based result elements.
#'
#' @param object A \code{microbTiSDA} object created by the package functions
#' @param max_preview Integer. Maximum number of rows and columns to preview for tabular results. Default is \code{5}
#' @param ... Additional arguments
#'
#' @export
#'
summary.microbTiSDA <- function(object, max_preview = 5, ...) {
  cat("Summary of MicrobTiSDA Object\n")
  cat("Class:", class(object), "\n")
  if (!is.null(object$params)) {
    cat("Parameters:\n")
    print(object$params)
  }
  if (!is.null(object$results)) {
    cat("Number of groups/items in results:", length(object$results), "\n")
    for (g in names(object$results)) {
      cat("\nGroup/Item:", g, "\n")
      if (is.data.frame(object$results[[g]])) {
        mat <- object$results[[g]]
        cat("Data frame dimensions:", dim(mat), "\n")
        cat("Preview (first", min(max_preview, nrow(mat)), "rows and columns):\n")
        print(mat[1:min(max_preview, nrow(mat)), 1:min(max_preview, ncol(mat))])
      } else if (is.list(object$results[[g]])) {
        for (nm in names(object$results[[g]])) {
          cat(" -", nm, ":\n")
          elem <- object$results[[g]][[nm]]
          if (is.matrix(elem) || is.data.frame(elem)) {
            cat("   Dimensions:", dim(elem), "\n")
            cat("   Preview:\n")
            print(elem[1:min(max_preview, nrow(elem)), 1:min(max_preview, ncol(elem))])
          } else {
            str(elem, max.level = 1)
          }
        }
      } else {
        str(object$results[[g]], max.level = 1)
      }
    }
  }
  invisible(object)
}


#' Summary Method for Design Objects
#'
#' @param object An object of class \code{Design}.
#' @param ... Additional arguments
#'
#' @export
summary.Design <- function(object, ...) {
  cat("Summary of Design object\n")
  cat("=================================\n")
  cat("Samples:", nrow(object$data), "\n")
  cat("Predictors:", ncol(object$data) - 2, "\n")
  cat("Time column: 'Time'\n")
  cat("ID column: 'ID'\n")
  if (!is.null(object$params$Group_var)) {
    cat("Grouping variable(s):", paste(object$params$Group_var, collapse = ", "), "\n")
    cat("Levels:", paste(colnames(object$data)[grepl("^interaction_terms", colnames(object$data))], collapse = ", "), "\n")
  }
  cat("---------------------------------\n")
  print(head(object$data, 5))
  invisible(object)
}


#' Summary Method for MicrobTiSDA Dynamic Interaction Visualization Objects
#'
#' Provides a concise summary of the dynamic interaction visualization results stored in a 'microbTiSDA_dynamic_vis' object.
#' The summary includes the number of nodes and edges in each group.
#'
#' @param object An object of class \code{microbTiSDA_dynamic_vis}, typically created by visualization functions in the \pkg{MicrobTiSDA} package.
#' @param ... Additional arguments.
#'
#' @export
summary.microbTiSDA_dynamic_vis <- function(object, ...) {
  cat("MicrobTiSDA Dynamic Interaction Visualization Summary\n")
  for (g in names(object$data)) {
    cat("Group:", g, "\n")
    cat("Number of nodes:", nrow(object$data[[g]]$nodes), "\n")
    cat("Number of edges:", nrow(object$data[[g]]$edges), "\n\n")
  }
  invisible(object)
}


#' Summary Method for MicrobTiSDA Spline Regression Objects
#'
#' Provides a concise summary for objects of class \code{MicrobTiSDA_spline_regression}. Prints the names of independent
#' variables and shows an example of the first fitted microbial feature.
#'
#' @param object An object of class \code{MicrobTiSDA_spline_regression}, usually created by spline regression functions in the \pkg{MicrobTiSDA} package.
#' @param ... Additional arguments.
#'
#' @export
summary.MicrobTiSDA_spline_regression <- function(object, ...) {
  cat("Summary of MicrobTiSDA Spline Regression Object\n")
  cat("Independent variables:", names(object$fitted_model), "\n")
  cat("Example of first fitted OTU:\n")
  print(object$fitted_model[[1]][[1]])
  invisible(object)
}


#' Summary Method for PredictedData Objects
#'
#' Provides a concise summary of a \code{PredictedData} object,
#' including number of groups, number of time points, number of features,
#' and the time range for each group.
#'
#' @param object A \code{PredictedData} object returned by \code{Pred.data}.
#' @param ... Additional arguments.
#'
#' @export
summary.PredictedData <- function(object, ...) {
  if (!inherits(object, "PredictedData")) {
    stop("Input must be a PredictedData object.")
  }

  summary_list <- lapply(object, function(df) {
    time_col <- "Predicted_Time"
    features <- setdiff(colnames(df), time_col)
    n_time <- nrow(df)
    n_features <- length(features)
    time_range <- range(df[[time_col]])
    data.frame(
      n_time_points = n_time,
      n_features = n_features,
      start_time = time_range[1],
      end_time = time_range[2]
    )
  })

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- names(object)

  cat("Summary of PredictedData object:\n")
  print(summary_df)

  invisible(summary_df)
}


#' Summary Method for MicrobTiSDA Cluster Objects
#'
#' Provides a summary of clustering results for objecs of class \code{MicrobTiSDA.cluster}. It reports the number of OTUs for each group, or
#' indicates if group was excluded.
#'
#' @param object An object of class \code{MicrobTiSDA.cluster}, typically produced by clustering functions in the \pkg{MicrobTiSDA} pacjage.
#' @param ... Additional arguments/
#' @export
summary.MicrobTiSDA.cluster <- function(object, ...) {
  cat("Summary of MicrobTiSDA Clustering Object\n")
  for (i in names(object$cluster_results)) {
    if (!is.null(object$cluster_results[[i]])) {
      cat(" Group:", i, " -> number of OTUs:",
          nrow(as.data.frame(object$predicted_data[[i]])), "\n")
    } else {
      cat(" Group:", i, " -> excluded (<3 OTUs)\n")
    }
  }
  invisible(object)
}

#' Summary Method for MicrobTiSDA ClusterCut OBjects
#'
#' Provides a summary of cluster cutting results for objects of class \code{MicrobTiSDA.clusterCut}. It reports the number of clusters for
#' each group, or indicates if no valid clustering was performend.
#'
#' @param object An object of class \code{MicrobTiSDA.clusterCut}, typically produced by cluster cutting function in the \pkg{MicrobTiSDA} package.
#' @param ... Additional arguments
#'
#' @export
summary.MicrobTiSDA.clusterCut <- function(object, ...) {
  cat("Summary of MicrobTiSDA Cluster Cutting\n")
  for (i in names(object$cut_assignments)) {
    if (!is.null(object$cut_assignments[[i]])) {
      cat(" Group:", i,
          " -> Clusters:", length(unique(object$cut_assignments[[i]]$cluster)), "\n")
    } else {
      cat(" Group:", i, " -> No valid clustering\n")
    }
  }
  invisible(object)
}

#' Summary Method for MicrobTiSDA Visual Objects
#'
#' Provides a summary of visualization results for objects of class \code{MicrobTiSDA.visual}. It reports, for each group, the number of cluster-speficif
#' visualization plot available.
#'
#' @param object An object of class \code{MicrobTiSDA.visual}, typically produced by visualization functions in the \pkg{MicrobTiSDA} package.
#' @param ... Additional arguments
#'
#' @export
summary.MicrobTiSDA.visual <- function(object, ...) {
  cat("Summary of MicrobTiSDA Visualizations\n")
  for (g in names(object$curves_plot)) {
    cat(" Group:", g,
        " -> Clusters:", length(object$curves_plot[[g]]), "\n")
  }
  invisible(object)
}

#' Summary method for RfBiomarker objects
#'
#' Provides a concise summary of the results from a Random Forest biomarker
#' analysis, including the number of selected OTUs, their names, prediction
#' results on training and test sets, and information about cross-validation.
#' @param object An object of class \code{"RfBiomarker"}.
#' @param ... Additional arguments
#'
#' @export
summary.RfBiomarker <- function(object, ...) {
  if (!inherits(object, "RfBiomarker")) stop("Input must be a RfBiomarker object.")

  cat("Summary of RfBiomarker object\n")
  cat("=============================\n\n")

  # 1. Number of selected OTUs
  num_otus <- nrow(object$OTU_importance)
  cat("Number of selected important OTUs:", num_otus, "\n\n")

  # 2. List selected OTUs
  cat("Selected OTUs:\n")
  print(rownames(object$OTU_importance))
  cat("\n")

  # 3. Training set predictions summary
  cat("Training set predictions summary:\n")
  if (!is.null(object$Predicted_results_on_train_set)) {
    print(table(object$Predicted_results_on_train_set))
  } else {
    cat("No training set predictions available.\n")
  }
  cat("\n")

  # 4. Test set predictions summary
  cat("Test set predictions summary:\n")
  if (!is.null(object$Predicted_results_on_test_set)) {
    print(table(object$Predicted_results_on_test_set))
  } else {
    cat("No test set predictions available.\n")
  }
  cat("\n")

  # 5. Cross-validation plot info
  cat("Cross-validation figure stored in the object.\n")
  cat("Use plot() to visualize the plot.\n")

  invisible(object)
}

#' Summary method for RegMESR objects
#'
#' Provides a concise summary of the results from a Regularized Multivariate
#' Exponential Spline Regression (RegMESR) analysis, including the number of
#' groups fitted, the total number of fitted models, and the model parameters.
#' @param object An object of class \code{"RegMESR"}.
#' @param ... Additional arguments
#'
#' @export
summary.RegMESR <- function(object, ...) {
  cat("Summary of RegMESR object\n")
  cat("---------------------------------\n")
  cat("Number of groups fitted: ", length(object$fitted_model), "\n")

  not_null_models <- sum(sapply(object$fitted_model, function(g) sum(!sapply(g, is.null))))
  cat("Total fitted models: ", not_null_models, "\n")

  cat("\nParameters:\n")
  print(object$parameters)

  invisible(object)
}


#' Summary method for PredictedDataMESR objects
#'
#' @param object An object of class PredictedDataMESR.
#' @param ... Additional arguments.
#'
#' @export
summary.PredictedDataMESR <- function(object, ...) {
  cat("Summary of PredictedDataMESR object\n")
  cat("===================================\n")
  cat("Number of groups:", length(object$groups), "\n\n")

  for (g in object$groups) {
    cat("Group:", g, "\n")
    pred_data <- object$data[[g]]
    cat("  Number of predictions:", nrow(pred_data), "\n")
    cat("  Features predicted:", paste(head(colnames(pred_data), 5), collapse = ", "))
    if (ncol(pred_data) > 5) cat(" ... (", ncol(pred_data), " total)\n", sep = "")
    cat("\n")
  }

  invisible(object)
}
