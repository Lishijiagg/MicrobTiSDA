# MicrobTiSDA: a flexible R package for inferring interspecies interactions and abundance dynamics in microbiome time series data
![Picture 1](https://github.com/user-attachments/assets/a16946af-fde8-4bad-8232-1585105d8dbb)


MicrobTiSDA is an R package specifically developed for the analysis of microbial community time-series data. It requires there types of input data:
1. Microbial compositional data table
2. Metadata
3. Microbial features annotations (e.g. OTU/ASV)

MicrobTiSDA consists of 7 functional modules, including:
- Input Module
- Data Preprocessing Module
- Species Intereaction Inference Module
- Regression Model Fitting Module
- Temporal Profiles CLustering Module
- Clustered Profiles Visualization Module
- Random Forest Classification Module

The core functionalities of MicrobTiSDA include:
1. Inferring internal species interactions within microbiomes based on discrete-time Lotka-Volterra model.
2. Constructing natural spline regression models to investigate temporal patterns of microbial abundances.
3. Identifying microbial features with similar (or opposite) temporal patterns through clustering analysis.

## Installation
Download the latest development code of MetaLonDA from GitHub using devtools
```r
library(devtools)
install_github("Lishijiagg/MicrobTiSDA")
```
