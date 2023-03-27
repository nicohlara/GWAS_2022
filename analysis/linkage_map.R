###Creates linkage maps for SunRILs populations
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-19

#set working directory
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

library(tidyverse)
library(gaston)
source("analysis/GWAS_functions.R")
library(ASMap)
library(onemap) #Mollinari collaborator mapping function

##Subset data for testing
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")
subset_geno <- select.inds(genotype, family=="UX1989")
#subset_geno <- select.snps(subset_geno, chr == "1A")

###ASMap testing
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


###onemap testing
body <- as.data.frame(t(as.matrix(subset_geno)))
#body <- body[1:1000,]

line_1 <- 'data type ri self'
line_2 <- c(dim(body)[2], dim(body)[1], 1, 1, 0)

body[body==0] <- 'a'
body[body==1] <-'-'
body[body==2] <-'b'
rownames(body) <- paste0('*',rownames(body))
#others <- data.frame('*', subset_geno@snps$id, )
sink('G:/My Drive/Nico_PhD/data/genotype/onemap_genotype/UX1989.raw')
cat(line_1, '\n')
cat(line_2, '\n')
cat(colnames(body), '\n')
sink()
#body2 <- cbind('A.H.B', body)
data.table::fwrite(cbind('A.B', body), 'G:/My Drive/Nico_PhD/data/genotype/onemap_genotype/UX1989.raw', sep=' ', quote = FALSE, row.names = TRUE, append = TRUE)
cat('*CHROM', ' ', subset_geno@snps$chr, '\n', append=TRUE, file='G:/My Drive/Nico_PhD/data/genotype/onemap_genotype/UX1989.raw')
cat('*POS', ' ', subset_geno@snps$pos, '\n', append=TRUE, file='G:/My Drive/Nico_PhD/data/genotype/onemap_genotype/UX1989.raw')

mapmaker <- read_onemap(dir = 'G:/My Drive/Nico_PhD/data/genotype/onemap_genotype', inputfile='UX1989.raw')
###also look into reading vcf in directly

##remove redundant markers to lower computational load
bins <- find_bins(mapmaker, exact=TRUE)
mapmaker_binned <- create_data_bins(mapmaker, bins)
##estimate recombination fractions
recomb_link <- rf_2pts(input.obj = mapmaker_binned, LOD = suggest_lod(mapmaker_binned))
##assign markers to linkage groups
recomb_mark_assig <- make_seq(recomb_link, 'all')
##these both do similar things: sort and condense into linkage groups
group_marker <- group(recomb_mark_assig)
chrom_group <- group_upgma(recomb_mark_assig, expected.groups = 21, inter = F)
##plot recombination in heatmap
LG1 <- make_seq(group_marker, 1)
LG1_rcd <- rcd(LG1, hmm=FALSE)
LG1_record <- record(LG1, hmm=FALSE)
LG1_ug <- ug(LG1, hmm=FALSE)

rf_graph_table(LG1_rcd)
