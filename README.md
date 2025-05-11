MicrobTiSDA: a flexible R package for inferring interspecies interactions and abundance dynamics in microbiome time series data
 
Author:   
Shijia Li

This package integrates the “Learning Interactions from Microbial Time Series” algorithm, which is based on the discrete-time Lotka-Volterra model, and supports the construction of natural spline regression models to analyze changes in microbial feature abundance over time. MicrobTiSDA not only enables the inference of species interaction networks and temporal patterns of species abundance but also provides methods to analyze the temporal dynamic similarity of species within microbial communities.
![Figure 1](https://github.com/user-attachments/assets/d4b93793-839b-4ef9-9ef8-85ee2dd6d1da)

## Installation
- Install the development version:
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

