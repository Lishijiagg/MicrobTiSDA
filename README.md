MicrobTiSDA: a flexible R package for inferring interspecies interactions and abundance dynamics in microbiome time series data
 
Author:   
Shijia Li

MicrobTiSDA integrates the “Learning Interactions from Microbial Time Series” algorithm, which is based on the discrete-time Lotka-Volterra model, and supports the construction of natural spline regression models to analyze changes in microbial feature abundance over time. MicrobTiSDA not only enables the inference of species interaction networks and temporal patterns of species abundance but also provides methods to analyze the temporal dynamic similarity of species within microbial communities.
![Figure 1](https://github.com/user-attachments/assets/d4b93793-839b-4ef9-9ef8-85ee2dd6d1da)

## Installation
- Install the development version of MicrobTiSDA:
```r
library(devtools)
install_github("Lishijiagg/MicrobTiSDA")
```

## Tutorial
Here is an example of applying MicrobTiSDA to an in vitro cultured aquatic microbiome dataset. The dataset was obtained from the study by [Fujita et al](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01474-5).

First, we should load MicrobTiSDA and related packages.
```r
library("MicrobTiSDA")
library("tidyr")
library("dplyr")
library("ggplot2")
```

In this dataset, the in vitro cultivation of aquatic microbiomes was conducted in eight replicate experiments. Each replicate involved continuous cultivation for 110 days, with daily sampling performed throughout the entire period. As a result, a total of 880 aquatic microbiome community samples were obtained. In total, 28 ASVs were detected across all samples.
```r
data("fujita.data")
data("fujita.meta")
data("fujita.taxa")
```

Next, we need to perform filtering on the feature table in this dataset. However, given the relatively small number of microbial features (ASVs) in this dataset, we chose not to filter out ASVs with low abundance or low prevalence.
```r
fujita_filt = Data.filter(Data = fujita.data,metadata = fujita.meta,OTU_counts_filter_value = 0,
                          OTU_filter_value = 0,Group_var = 'replicate.id')
```

Typically, the filtered data need to be transformed prior to analysis. MicrobTiSDA provides a modified centered log-ratio (MCLR) transformation for this purpose. This transformation is performed separately for each sampled individual or experimental group; therefore, users are required to specify the grouping variable from the metadata. In this dataset, the variable representing the eight replicate experiments is “replicate.id”.
```r
fujita_trans = Data.trans(Data = fujita_filt,metadata = fujita.meta,Group_var = 'replicate.id')
```

After filtering and transforming the dataset, MicrobTiSDA can be used to infer species interactions within each individual subject or replicate microbiome. By integrating the ["Learning Interactions from Microbial Time-Series" (LIMITS) framework](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102451), MicrobTiSDA is able to infer sparse interaction networks from microbiome time-series data. This step may take a considerable amount of time to complete, depending on the size of the dataset.
```r
fujita_interact = Spec.interact(Data = fujita_trans,metadata = fujita.meta,Group_var = 'replicate.id',num_iterations = 10)
```

Next, we can visualize the results of the species interaction inference. The “Interact.vis” function generates network plots of species interactions for each replicate aquatic microbiome experiment.
```r
fujita_interact_vis = Interact.vis(Interact_data = fujita_interact,count_data = fujita_trans,metadata = fujita.meta,
                                   Subject_ID = 'replicate.id',legend_title_size = 12,legend_text_size = 10,
                                   label_distance = 0.35,label_font_size = 4,title_size = 15)
```

```r
fujita_interact_vis$Rep.1
```
![image](https://github.com/user-attachments/assets/83c528e7-8e14-4cd7-a954-68eb515007b9)

After completing species interaction inference for the aquatic microbiome, we can proceed to analyze the temporal abundance patterns of microbial features using MicrobTiSDA. To support this analysis, MicrobTiSDA provides a method for fitting regression models to the abundance trajectories of microbial features. To begin, we first need to construct a design matrix for modeling temporal trends for each microbial feature. This can be accomplished using the Design function, where users must specify the variables in the metadata that represent the group ID for each subject or replicate, the sample ID, and the time point of each sample.
```r
fujita_design = Design(metadata = fujita.meta,Group_var = 'replicate.id',Sample_ID = 'timeChar',Pre_processed_Data = fujita_trans,Sample_Time = 'time')
```

When fitting regression models for each microbial feature, MicrobTiSDA employs natural spline regression. Users can manually specify the positions of spline knots by providing a vector—for example, c(10, 20, 30) sets the 10th, 20th, and 30th days in the time series as knot positions. Alternatively, MicrobTiSDA can automatically determine the optimal number and placement of knots using generalized cross-validation, based on a user-defined maximum number of knots.
```r
fujita_model = Reg.SPLR(Data_for_Reg = fujita_design,
                      pre_processed_data = fujita_trans,
                      max_Knots = 5,
                      unique_values = 5)
```
here, it is important to account for the sparsity and zero-inflated nature of microbiome data. The "Reg.SPLR" function includes the "unique_values" parameter, which specifies the minimum number of non-zero abundance values required across time-series samples for a microbial feature to be included in regression modeling. By default, this threshold is set to at least 5 non-zero values.

Based on the fitted regression models, we can predict the temporal abundance patterns of each microbial feature. These features can then be clustered according to the similarity of their abundance temporal patterns, using correlation distance as the clustering metric.
```r
fujita_model_pred = Pred.data(Fitted_models = fujita_model,
                            metadata = fujita.meta,
                            Group = 'replicate.id',
                            Sample_Time = 'time',
                            time_step = 1)

fujita_model_clust = Data.cluster(predicted_data = fujita_model_pred,
                                clust_method = 'average',
                                dend_title_size = 12,
                                font_size = 3)
```
![image](https://github.com/user-attachments/assets/612945b3-970c-4dab-9077-07db85febbcc)

We can also identify the optimal clusters of microbial features based on the clustering results. Here we select features with correlation distance less than 0.15 for clustering.
```r
fujita_clust_results = Data.cluster.cut(cluster_outputs = fujita_model_clust,cut_height = 0.15,font_size = 4)

fujita_clust_results$cluster_figures$Rep.1
```
![image](https://github.com/user-attachments/assets/bf8cb876-94e9-42f5-8808-df3f67d8a10b)

Finally, we can visualize the abundance patterns of microbial features within each selected cluster.
```r
fujita_model_vis = Data.visual(cluster_results = fujita_clust_results,
                             cutree_by = 'height',
                             cluster_height = rep(0.15,8),
                             predicted_data = fujita_model_pred,
                             Design_data = fujita_design,
                             pre_processed_data = fujita_trans,
                             plot_dots = TRUE,
                             figure_x_scale = 20,
                             Taxa = fujita.taxa,
                             plot_lm = FALSE,
                             legend_title_size = 10,
                             legend_text_size = 10,axis_x_size = 10,axis_title_size = 10,axis_y_size = 10)

fujita_model_vis$Rep.1[[4]]
```
![image](https://github.com/user-attachments/assets/3c311587-1ddf-48bd-97e3-64f38ca37b0a)
