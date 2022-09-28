###Creates linkage maps for SunRILs populations
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-19

#set working directory
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

library(tidyverse)
library(gaston)
source("analysis/GWAS_functions.R")
library(ASMap)

genotype <- read.bed.matrix("data/SunRILs_filtered_2022")
subset_geno <- select.inds(genotype, family=="UX1989")
#subset_geno <- select.snps(subset_geno, chr == "1A")
geno <- as.data.frame(t(as.matrix(subset_geno)), stringAsFactors = FALSE)
geno[geno==0] <- "A"
geno[geno==1] <- "a"
geno[geno==2] <- "B"
test_map <- mstmap(geno, bychr=TRUE, "RIL6")
heatMap(test_map)

map_est <- quickEst(test_map)
map_est <- subset(map_est, chr = names(nmar(map_est))[6:15])
plot.map(test_map)


map2 <- quickEst(mapDH, map.function = "kosambi")
map2 <- subset(map2, chr = names(nmar(map2))[6:15])
plot.map(map1, map2)


#visualization test
map_test <- subset(geno)
alignCross(geno, chr = "1A", maps = row.names(map_test) )

pull.map(test_map)

#data(mapDH, package = "ASMap")
