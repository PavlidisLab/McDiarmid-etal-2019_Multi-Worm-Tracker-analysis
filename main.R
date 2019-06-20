#!/usr/bin/env Rscript

print("Generating data frame")
source("preprocessing_scripts/create_dataframe.R")

print("Processing phenotypes")
source("preprocessing_scripts/process.R")

print("Computing statistics")
source("preprocessing_scripts/t_test.R")

print("Clustering phenotypes")
source("preprocessing_scripts/clustering-t_stat.R")

print("Generating t-SNE")
source("preprocessing_scripts/tsne.R")
