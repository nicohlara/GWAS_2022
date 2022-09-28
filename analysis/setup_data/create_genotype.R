###Create and clean genotype files for further analysis
##Created by Nicolas A. H. Lara
##Last edit: 2022-9-12

###SET WORKING DIRECTORY
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
library(gaston)

###READ IN IMPUTED VCF FILE, CLEAN UP DATA
##IMPUTED USING BEAGLE5 ###only has some populations
#genotype <- read.vcf("data/used_vcf/SunRILs_2021_postimp_filt.vcf.gz", convert.chr = F)
##IMPUTED USING NOAH'S ALGORITHM###
genotype <- read.vcf("data/SunRILs_combGenos_fakeHeader_filtered_imp.vcf.gz", convert.chr = F)
#filter out parents and thin on LD
#genotype <- select.inds(genotype, grepl("^UX", id))
genotype <- LD.thin(genotype, threshold = .8, max.dist = 350e6) #made for humans, Noah suggested changing to ~350MB 
genotype_matrix <- as.matrix(genotype)
#reformat ids to get rid of extra text
genotype@ped$id <- gsub("-NWG", "", genotype@ped$id)
genotype@ped$id <- gsub("-A+", "", genotype@ped$id)
genotype@ped$id <- gsub("-NEG", "", genotype@ped$id)
#filter for resequenced and duplicated lines
genotype <- select.inds(genotype, !grepl("-A-", id))
genotype <- select.inds(genotype, !duplicated(id))
##Add family marker for selection

a <- data.frame(ID = genotype@ped$famid) %>% separate(ID, c("family"), extra = "drop")
parents <- c("GA00190", "MPV57", "GA001138", "SCTX98", "NC08", "HILLIARD", "GA05450", "SS8641", "NC8248", "LA09264C", "GA06493", "AGS2000", "TX12D4896", "LA95135", "LA0964C", "ARGA051160") 
for (i in c(1:length(a[,1]))) {
	if (a[i,1] %in% parents) {
		a[i,1] <- "Parents" 
	}
}
genotype@ped$family <- a[,1]
# ###filter by planted entries
# genotype_sync <- select.inds(genotype, id %in% phenotype$Entry)
# ###save resulting genotype file
# write.bed.matrix(genotype_sync, "data/SunRILs_filtered_2022")

write.bed.matrix(genotype, "data/SunRILs_filtered_2022")

