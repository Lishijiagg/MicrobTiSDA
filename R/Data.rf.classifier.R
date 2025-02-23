#' @title Random Forest classification and Biomarker selection for OTU/ASV data
#'
#' @description
#' This function implements a random forest classification model tailored for OTU/ASV datasets.
#'     It performs data filtering, model training, performance evaluation, cross-validation, and
#'     biomarker (important microbial features) selection based on Mean Decrease Accuracy.
#'
#' @details
#' The function processes the input OTU count data and corresponding metadata in several steps:
#' \enumerate{
#'   \item \strong{Data Filtering and Preparation:} If a minimum count threshold (`OTU_counts_filter_value`) is provided,
#'   OTUs with total counts below this value are removed. The OTU table is then transposed and merged with the metadata,
#'   where a specific column (specified by `Group`) indicates the group labels.
#'   \item \strong{Data Partitioning:} The combined dataset is split into training and testing subsets based on the proportion
#'   specified by `train_p`.
#'   \item \strong{Model Training:} A random forest classifier is trained on the training data. The function computes the
#'   margin scores for the training samples, which are plotted to visualize the model’s confidence.
#'   \item \strong{Performance Evaluation:} Predictions are made on both training and testing datasets. Confusion matrices
#'   are generated to compare the actual versus predicted classes.
#'   \item \strong{Feature Importance and Cross-Validation:} OTU importance is assessed using Mean Decrease Accuracy.
#'   Repeated k-fold cross-validation (default 10-fold repeated `reps` times) is performed to determine the optimal number
#'   of OTUs (biomarkers). A cross-validation error curve is plotted, and the user is prompted to input the best number
#'   of OTUs based on the plot.
#'   \item \strong{Biomarker Selection:} Based on the user input, the top OTUs are selected as potential biomarkers, and
#'   the corresponding OTU table is returned along with all the diagnostic outputs.
#' }
#'
#' @param raw_data A numeric matrix or data frame of counts data with OTUs/ASVs as rows and samples as columns.
#' @param metadata A data frame. Containing information about all samples, including at least the grouping of all samples as well as
#'     individual information (\code{Group} and \code{ID}), the sampling \code{Time} point for each sample, and other relevant information.
#' @param train_p A positive decimal. Indicating the percentage of data that goes to training. For example, when
#'     \code{train_p = 0.7}, 70% samples were randomly selected as training dataset. More information see
#'     \code{\link[randomForest]{rfcv}}.
#' @param Group A string that specifies the columns in the \code{metadata} for grouping the temporal series samples.
#' @param OTU_counts_filter_value An integer, indicating the sum of the minimum abundances of OTUs/ASVs in all samples. If the sum
#'     of the abundances that OTU/ASV is below the given positive integer threshold, the OTU/ASV is excluded, and vice versa, it
#'     is retained. The default is \code{NA}.
#' @param reps An integer. The number of replications for cross-validation. By default, \code{reps = 5}. More details see
#'     \code{\link[randomForest]{rfcv}}.
#' @param cv_fold An integer. Number of folds in the cross-validation. By default, \code{cv_fold = 10}. see \code{\link[randomForest]{rfcv}}
#' @param title_size Numeric value for the font size of plot titles. Defaults to 10.
#' @param axis_title_size Numeric value for the font size of axis titles. Defaults to 8.
#' @param legend_title_size Numeric value for the font size of legend titles. Defaults to 8.
#' @param legend_text_size Numeric value for the font size of legend text. Defaults to 6.
#'
#' @return A list containing the following elements: (1) A data frame of the selected OTU/ASV table containing the top biomarkers.
#'     (2) Predicted class labels for the training dataset. (3) Predicted class labels for the testing dataset. (4) A confusion
#'     matrix comparing actural versus predicted classes in the training set. (5) A confusion matrix comparing actural versus predicted
#'     classes in the testing set. (6) A ggplot object displaying the margin scores for the training dataset. (7) A data frame of OTU
#'     importance metrics from the random forest model, ranked by Mean Decrease Accuracy. (8) A duplicate of the testing set confusion
#'     matrix (for compatibility). (9) A ggplot object of the cross-validation error curve with a vertical line indicating the user-selected
#'     optimal number of OTUs/ASVs.
#'
#' @importFrom caret createDataPartition
#' @importFrom randomForest randomForest
#' @importFrom randomForest margin
#' @importFrom randomForest rfcv
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 melt
#' @importFrom splines ns
#' @export
#' @author Shijia Li
#'
#' @examples
#' \dontrun{
#' # Example OTU count data (20 OTUs x 10 samples)
#' set.seed(123)
#' otu_data <- matrix(sample(0:100, 200, replace = TRUE), nrow = 20)
#' colnames(otu_data) <- paste0("Sample", 1:10)
#' rownames(otu_data) <- paste0("OTU", 1:20)
#'
#' # Example metadata with group labels
#' metadata <- data.frame(Group = rep(c("Control", "Treatment"), each = 5))
#'
#' # Run the classifier
#' result <- Data.rf.classifier(raw_data = otu_data,
#'                              metadata = metadata,
#'                              train_p = 0.7,
#'                              Group = "Group",
#'                              OTU_counts_filter_value = 50)
#' }
Data.rf.classifier = function(raw_data,
                          metadata,
                          train_p,
                          Group,
                          OTU_counts_filter_value = NA,
                          reps = 5,
                          cv_fold = 10,
                          title_size = 10,
                          axis_title_size = 8,
                          legend_title_size = 8,
                          legend_text_size = 6) {

  otu = raw_data

  if (!is.na(OTU_counts_filter_value)) {
    otu = otu[which(rowSums(otu) >= OTU_counts_filter_value),]
  } else {
    otu = otu
  }
  otu = as.data.frame(t(otu))
  meta = metadata

  otu$group = meta[,Group]
  otu$group = as.factor(otu$group)

  # devide otu table into train dataset and test dataset
  set.seed(123)
  select_train = caret::createDataPartition(otu$group,p = train_p,list = FALSE)
  otu_train = otu[select_train,]
  otu_train$group = as.factor(otu_train$group)
  otu_test = otu[-select_train,]

  # fitting a classification model on the training dataset using randomForest
  set.seed(123)
  otu_train.forest = randomForest::randomForest(group ~ ., data = otu_train, importance = TRUE)
  print('==================================================================================')
  print(otu_train.forest)
  print('==================================================================================')
  margin_scores <- randomForest::margin(otu_train.forest, otu_train$group)
  margin_df <- data.frame(
    Group = names(margin_scores),
    Margin_Score = as.vector(margin_scores)
  ) %>% arrange(Margin_Score) %>% mutate(Sample = row_number())
  oob_error = round(otu_train.forest$err.rate[nrow(otu_train.forest$err.rate),"OOB"],3)
  p_margin = ggplot(margin_df, aes(x = Sample, y = Margin_Score,color = Group)) +
    geom_point() +
    labs(title = 'Margin Scores in Train-Dataset',
         x = 'Number of Samples',
         y = 'Margin Score') +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          title = element_text(size = title_size),
          axis.title = element_text(size = axis_title_size),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size))+
    annotate("text",x = Inf,y = -Inf,hjust = 1.1,vjust = -0.5,
             label = paste("Out of Bag error = ",oob_error,sep = ""),
             size = 4,color = "black")

  ## classifier performance testing

  # training set self-testing
  train_predict = predict(otu_train.forest, otu_train)
  compare_train = table(train_predict,otu_train$group, dnn = c('Actural','Predicted'))

  # classifier evaluation by testing dataset
  test_predict = predict(otu_train.forest,otu_test)
  compare_test = table(otu_test$group, test_predict, dnn = c('Actural','Predicted'))

  # biomarker selections
  importance_otu = otu_train.forest$importance
  importance_otu = as.data.frame(importance_otu)
  # rank important OTUs based on Mean Decrease Accuracy
  importance_otu = importance_otu[order(importance_otu$MeanDecreaseAccuracy,
                                        decreasing = TRUE),]

  ## cross-validation for selecting the best number of biomarkers
  # repeat 10-fold cross-validations for 5 times
  set.seed(123)
  otu_train.cv = replicate(reps, randomForest::rfcv(otu_train[-ncol(otu_train)],
                                                    otu_train$group,cv.fold = cv_fold,step = 1.5),
                           simplify = FALSE)

  # extract validation results for plotting
  otu_train.cv = data.frame(sapply(otu_train.cv,'[[','error.cv'))
  otu_train.cv$otus = rownames(otu_train.cv)
  otu_train.cv = reshape2::melt(otu_train.cv,id = 'otus')
  otu_train.cv$otus = as.numeric(as.character(otu_train.cv$otus))

  p = ggplot(otu_train.cv,aes(otus,value)) +
    geom_smooth(se = FALSE, method = 'glm', formula = y~splines::ns(x,6))+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          title = element_text(size = title_size),
          axis.title = element_text(size = axis_title_size),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size)) +
    labs(title = 'Cross-validation curve',
         x = 'Number of Features',
         y = 'Cross-validation error')
  print(p)

  print('=========Classification results of the classifier on the training dataset=========')
  print(compare_train)
  print('==================================================================================')
  print('==========Classification results of the classifier on the testing dataset=========')
  print(compare_test)
  print('==================================================================================')

  features_number = readline(paste(
    "According to the cross-validation results, the best selected number of OTUs could be: "))
  p = p + geom_vline(xintercept = as.numeric(features_number),linetype = 'dashed')+
    scale_x_continuous(breaks = c(as.numeric(features_number)))
  print(p)

  otu_select = rownames(importance_otu)[1:features_number]
  important_otu_table = otu[otu_select]
  important_otu_table = as.data.frame(t(important_otu_table))

  output_result = list(important_otu_table,
                       train_predict,
                       test_predict,
                       compare_train,
                       compare_test,
                       p_margin,
                       importance_otu,
                       compare_test,
                       p)
  names(output_result) = c('Important_OTU_table',
                           'Predicted_results_on_train_set',
                           'Predicted_results_on_test_set',
                           'Traindata_confusion_matrix',
                           'Testdata_confusion_matrix',
                           'Margin_scores_train',
                           'OTU_importance',
                           'Confusion_matrix_testset',
                           'cross-validation')

  return(output_result)
}
