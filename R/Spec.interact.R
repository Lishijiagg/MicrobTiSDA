#' @title Species Interaction Inferrences
#' @description
#' This function describes interspecies interactions based on the discrete-time Lotka-Volterral model.
#' @details
#' This function implements the discrete-time Lotka-Volterra model to characterize species interactions in microbiome time-series data.
#'     The model describes the abundance (MCLR transformed) \eqn{x_{ni}} of species \eqn{i} for subject \eqn{n} at time \eqn{t+\Delta t} as:
#'     \deqn{x_{ni} (t+\Delta t) = \eta_{ni} (t) x_{ni} (t) \exp\left(\Delta t \sum_j c_{nij} (x_{nj} (t) - <x_{nj}>) \right)}
#'     where \eqn{<x_{nj}>} represents the equilibrium abundance of species \eqn{j}, typically defined as the
#'     median abundance across samples from the same subject; \eqn{c_{nij}} denotes the interaction coefficient of species
#'     \eqn{j} on species \eqn{i}; and \eqn{\eta_{ni} (t)} accounts for log-normally distributed stochastic effects. For
#'     computational simplicity, stochastic effects are ignored, \eqn{\Delta t} is set to 1. Taking the natural logarithm yealds:
#'     \deqn{\ln x_{ni} (t+1) - \ln x_{ni} (t) = \sum_j c_{nij} (x_{nj} (t) - <x_{nj}>)}
#'     To improve sparsity and interpretability, the LIMITS algorithm is applied, incorporating stepwise regression and bagging.
#'     First, 50% of the samples are randomly selected as the training set while the rest serve as the test set. An initial regression
#'     model includes only the self-interaction term:
#'     \deqn{\ln x_{ni} (t+1) - \ln x_{ni} (t) = c_{nii} (x_{ni} (t) - <x_{ni}>)}
#'     Stepwise regression then iteratively adds species interaction terms from a candidate set \eqn{S}, forming:
#'     \deqn{\ln x_{ni} (t+1) - \ln x_{ni} (t) = c_{nii} (x_{ni} (t) - <x_{ni}>) + \sum_{j \in S} c_{nij} (x_{nj} (t) - <x_{nj}>)}
#'     The inclusion of a new term is determined based on the improvement in mean squared error (MSE) on the test set:
#'     \deqn{\theta = \frac{\text{MSE}_{\text{before}} - \text{MSE}_{\text{after}}}{\text{MSE}_{\text{before}}}}
#'     If \eqn{\theta} exceeds a predefined threshold (default \( 10^{-3} \)), the species is included. Bagging is performed over \eqn{B}
#'     iterations by repeating the random splitting and stepwise regression, enhancing robustness. The final interaction coefficient matrix is computed as:
#'     \deqn{c_{nij} = \text{median}(c_{nij}^{(1)}, c_{nij}^{(2)}, ..., c_{nij}^{(B)})}
#'     This approach refines the inferred species interactions while ensuring sparsity.
#'
#'
#' @param Data A matrix or data frame of the transformed species abundance data.
#' @param metadata A data frame. Containing information about all samples, including at least the grouping of all samples as well as
#'     individual information (\code{Group} and \code{ID}), the sampling \code{Time} point for each sample, and other relevant information.
#' @param Group_var A character string specifying the column name in \code{metadata} that defines the groups for analysis.
#' @param abund_centered_method A character string indicating the method to compute species equilibrium abundance.
#'     Accepted values are \code{median} (default) and \code{mean}.
#' @param num_iterations An integer specifying the number of bagging iterations for the iterative variable selection process. Default is 10.
#' @param error_threshold A numeric value representing the relative error improvement threshold for adding new predictors during bagging iteration.
#'     Default is 1e-3.
#' @param pre_error A numeric value specifying the initial (large) error used for comparison in the iterative procedure. Default is 10000.
#'
#' @return A list with an element for each group defined by \code{Group_var}. Each element is a list containing:
#' \describe{
#'   \item{interaction_matrices}{A three-dimensional array of estimated interaction coefficients with dimensions corresponding to features
#'   \eqn{\times} features \eqn{\times} iterations.}
#'   \item{final_interaction_matrix}{A two-dimensional matrix of interaction coefficients obtained by taking the median over the iterations.}
#' }
#' @export
#' @author Shijia Li
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assume 'abundance_data' is a matrix of species abundances (samples x species)
#' # and 'sample_metadata' is a data frame with sample information including a grouping variable "Treatment".
#' results <- Spec.interact(Data = abundance_data,
#'                          metadata = sample_metadata,
#'                          Group_var = "Treatment",
#'                          abund_centered_method = "median",
#'                          num_iterations = 20,
#'                          error_threshold = 1e-3,
#'                          pre_error = 10000)
#' }
Spec.interact <- function(Data, metadata, Group_var, abund_centered_method = 'median',
                           num_iterations = 10, error_threshold = 1e-3, pre_error = 10000) {
  results <- list()
  Data = as.data.frame(t(Data))

  for (g in unique(metadata[[Group_var]])) {
    group_info <- metadata[metadata[[Group_var]] == g, ]
    otu_table <- as.data.frame(t(Data[, rownames(group_info)] + 1))

    equilibrium_abundance <- if (abund_centered_method == "median") {
      apply(otu_table, 2, median)
    } else if (abund_centered_method == "mean") {
      apply(otu_table, 2, mean)
    } else {
      stop("abund_centered_method must be 'median' or 'mean'")
    }

    design_matrix_list <- list()
    data_vector_list <- list()
    for (i in seq_along(otu_table)) {
      log_diff_vector <- diff(log(otu_table[[i]]))
      design_matrix <- as.matrix(otu_table[-nrow(otu_table), ] - equilibrium_abundance)
      design_matrix_list[[i]] <- design_matrix
      data_vector_list[[i]] <- log_diff_vector
    }

    interaction_matrices <- array(0, dim = c(ncol(otu_table), ncol(otu_table), num_iterations))

    for (iter in 1:num_iterations) {
      set.seed(iter)
      sample_indices <- sample(1:(nrow(otu_table) - 1), size = floor((nrow(otu_table) - 1) / 2))
      training_data <- otu_table[sample_indices, ]
      test_data <- otu_table[-sample_indices, ]

      for (i in seq_along(otu_table)) {
        active_set <- c(i)
        inactive_set <- setdiff(seq_len(ncol(otu_table)), i)
        data_vector <- data_vector_list[[i]][sample_indices]
        test_data_vector <- data_vector_list[[i]][-sample_indices]

        prev_error <- pre_error
        while (length(inactive_set) > 0) {
          candidate_errors <- sapply(inactive_set, function(j) {
            candidate_set <- c(active_set, j)
            design_matrix <- as.matrix(training_data[, candidate_set] - equilibrium_abundance[candidate_set])
            if (nrow(design_matrix) < ncol(design_matrix)) return(Inf)

            candidate_coefs <- MASS::ginv(design_matrix) %*% data_vector
            if (any(is.na(candidate_coefs) | is.infinite(candidate_coefs))) return(Inf)

            pred_values <- as.matrix(test_data[, candidate_set] - equilibrium_abundance[candidate_set]) %*% candidate_coefs
            error <- suppressWarnings(mean((pred_values - test_data_vector)^2))

            if (is.na(error) || is.infinite(error)) return(Inf)
            return(error)
          })

          if (length(candidate_errors) == 0 || all(is.infinite(candidate_errors))) break
          best_j <- inactive_set[which.min(candidate_errors)]
          best_error <- min(candidate_errors)

          if (!is.na(prev_error) && !is.na(best_error) && prev_error > 0 && best_error >= 0) {
            if ((prev_error - best_error) / prev_error > error_threshold) {
              active_set <- c(active_set, best_j)
              inactive_set <- setdiff(inactive_set, best_j)
              prev_error <- best_error
            } else {
              break
            }
          } else {
            break
          }
        }

        final_design_matrix <- as.matrix(training_data[, active_set] - equilibrium_abundance[active_set])
        final_data_vector <- data_vector_list[[i]][sample_indices]
        final_coefs <- MASS::ginv(final_design_matrix) %*% final_data_vector
        interaction_matrices[active_set, i, iter] <- final_coefs
      }
    }

    final_interaction_matrix_median <- apply(interaction_matrices, c(1, 2), median, na.rm = TRUE)
    colnames(final_interaction_matrix_median) <- colnames(otu_table)
    rownames(final_interaction_matrix_median) <- colnames(otu_table)

    results[[g]] <- list(interaction_matrices = interaction_matrices, final_interaction_matrix = final_interaction_matrix_median)
  }

  return(results)
}
