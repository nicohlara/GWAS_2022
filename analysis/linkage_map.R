###Creates linkage maps for SunRILs populations
#Created by Nicolas A. H. Lara
#Last edit: 2022-8-24

#set working directory
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

library(tidyverse)
library(gaston)
source("analysis/GWAS_functions.R")

genotype <- read.bed.matrix("data/SunRILs_filtered_2022")