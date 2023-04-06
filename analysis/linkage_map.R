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
#genotype <- read.bed.matrix("data/SunRILs_filtered_2022")
#subset_geno <- select.inds(genotype, family=="UX1989")
#subset_geno <- select.snps(subset_geno, chr == "1A") 

genotype <- read.vcf('G:/My Drive/Nico_PhD/data/genotype/SunRILs_2021_UX1992_RILs.vcf', convert.chr = FALSE)
a <- data.frame(ID = genotype@ped$famid) %>% separate(ID, c("family"), extra = "drop")
parents <- c("GA00190", "MPV57", "GA001138", "SCTX98", "NC08", "HILLIARD", "GA05450", "SS8641", "NC8248", "LA09264C", "GA06493", "AGS2000", "TX12D4896", "LA95135", "LA0964C", "ARGA051160") 
for (i in c(1:length(a[,1]))) {
  if (a[i,1] %in% parents) {
    a[i,1] <- "Parents" 
  }
}
genotype@ped$family <- a[,1]

subset_geno <- select.inds(genotype, family == 'UX1992')
subset_geno <- select.snps(subset_geno, callrate >= .9 & maf >=.3 & chr != 'UN')


#############################ASMap testing########################################
###testing subset
subset_geno <- select.snps(subset_geno, chr == '6A')

###restructuring data to qtl object
geno <- as.data.frame(t(as.matrix(subset_geno)), stringAsFactors = FALSE)
geno[geno==0] <- "A"
geno[geno==1] <- "X"
geno[geno==2] <- "B"
geno[is.na(geno)] <- "U"

###creating basic object for analysis and preprocessing
RIL_map <- mstmap(geno, bychr=TRUE, "RIL6")
RIL_map_f <- pullCross(RIL_map, type="co.located")
#RIL_map_f <- pullCross(RIL_map_f, type="seg.distortion")

###remake map
RIL_map_f <- mstmap(RIL_map_f, bychr = FALSE)



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

# ###onemap testing
# body <- as.data.frame(t(as.matrix(subset_geno)))
# #body <- body[1:1000,]
# 
# line_1 <- 'data type ri self'
# line_2 <- c(dim(body)[2], dim(body)[1], 1, 1, 0)
# 
# line_name <- 'UX1992'
# 
# file_name <- paste('G:/My Drive/Nico_PhD/data/genotype/onemap_genotype/', line_name, '.raw', sep="")
# body[body==0] <- 'a'
# body[body==1] <-'-'
# body[body==2] <-'b'
# rownames(body) <- paste0('*',rownames(body))
# #others <- data.frame('*', subset_geno@snps$id, )
# sink(file_name)
# cat(line_1, '\n')
# cat(line_2, '\n')
# cat(colnames(body), '\n')
# sink()
# #body2 <- cbind('A.H.B', body)
# data.table::fwrite(cbind('A.B', body), file_name, sep=' ', quote = FALSE, row.names = TRUE, append = TRUE)
# cat('*CHROM', ' ', subset_geno@snps$chr, '\n', append=TRUE, file=file_name)
# cat('*POS', ' ', subset_geno@snps$pos, '\n', append=TRUE, file=file_name)
# 
# mapmaker <- read_onemap(dir = 'G:/My Drive/Nico_PhD/data/genotype/onemap_genotype', inputfile=paste(line_name, '.raw', sep=""))
# ###also look into reading vcf in directly

vcf_mapmaker <- onemap_read_vcfR(vcf='G:/My Drive/Nico_PhD/data/genotype/SunRILs_2021_UX1992_RILs.vcf', parent1='HILLIARD', parent2='GA06493-13LE6', cross='ri self')
mapmaker <- filter_missing(vcf_mapmaker, threshold = .10, by='markers')
##remove redundant markers to lower computational load
bins <- find_bins(mapmaker, exact=TRUE)
mapmaker_binned <- create_data_bins(mapmaker, bins)
##estimate recombination fractions
recomb_link <- rf_2pts(input.obj = mapmaker_binned, LOD = suggest_lod(mapmaker_binned))
for (chrom in unique(mapmaker_binned$CHROM)) {
  chrom <- '1A'
  ##assign markers to linkage groups
  recomb_mark_assig <- make_seq(recomb_link, chrom)#'all')
  ##these both do similar things: sort and condense into linkage groups
  group_marker <- group(recomb_mark_assig, LOD = suggest_lod(recomb_mark_assig))
  #chrom_group <- group_upgma(recomb_mark_assig, inter = T, comp.mat = TRUE)
  ##plot recombination in heatmap
  LG1 <- make_seq(group_marker, 1)
  #LG1_rcd <- rcd(LG1, hmm=FALSE)
  #LG1_record <- record(LG1, hmm=FALSE)
  LG1_ug <- ug(LG1, hmm=FALSE)

  rf_graph_table(LG1_ug, mrk.axis='names')
}