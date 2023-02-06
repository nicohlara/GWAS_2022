#####Using the qtl2 package, create a cross2 object and perform analysis on it
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-13

###SET WORKING DIRECTORY
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
library(tidyverse)
library(gaston)
library(qtl2)
library(qtl2convert)

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")

###creating cross2 files
# genotype_subset <- select.inds(genotype, id %in% phenotype$Entry)
# geno <- as.matrix(genotype_subset)
# #geno <- t(geno
# allele_codes <- as.matrix(genotype@snps[c('A1', 'A2')])
# rownames(allele_codes) <- genotype@snps$id
# for (i in colnames(geno)) {
#   geno[,i][geno[,i]==0] <- allele_codes[i,1]
#   geno[,i][geno[,i]==1] <- paste(allele_codes[i,1], allele_codes[i,2], sep="")
#   geno[,i][geno[,i]==2] <- allele_codes[i,2]
# }
# output_codes <- c('-', "A", "H", "B")
# geno <- t(encode_geno(t(geno), allele_codes))
# for (i in colnames(geno)) {
#   if (length(unique(geno[,i])) <= 1) {geno <- geno[,colnames(geno)!=i]}
# }
#write.csv(geno, "output/data/SunCross/SunRILs_geno.csv")
# gmap <- subset(genotype@snps, select=c("id", "chr", "pos"))
# names(gmap)[names(gmap)=="id"] <- "marker"
# write.csv(gmap, "output/data/SunCross/SunRILs_gmap.csv", row.names=FALSE)
pheno <- subset(phenotype, Location == "Raleigh", select= -c(Location, Cross_ID) )#select=c("Entry", "Location", "Height"))
# #pheno_K <- subset(pheno, Location == "Kinston", select=c("Entry", "Height"))
# #pheno_K <- rename(pheno_K, Kinston = Height)
# #pheno_R <- subset(pheno, Location == "Raleigh", select=c("Entry", "Height"))
# #pheno_R <- rename(pheno_R, Raleigh = Height)
# #pheno <- merge(pheno_R, pheno_K, by="Entry")
write.csv(pheno, "output/data/SunCross/SunRILs_pheno.csv", row.names=FALSE)
# write.csv(phenotype, "output/data/SunCross/SunRILs_pheno.csv", row.names=FALSE)
# #phenocovar <- data.frame(pheno=c("Kinston", "Raleigh"), location = c("Kinston", "Raleigh"))
# #write.csv(phenocovar, "/Users/nico/Desktop/cross2_file_test/SunRILs_phenocovar.csv", row.names=FALSE)
# covar <- distinct(subset(phenotype, select=c("Entry", "Cross_ID")))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995")))
# write.csv(covar, "output/data/SunCross/SunRILs_covar.csv", row.names=FALSE)
#cross_parents <- read.delim('output/data/SunCross/SunRILs_cross_parents.csv', sep=",")

#write.csv(cross_info, "output/data/SunCross/SunRILs_cross_info.csv")
SunCross <- read_cross2("output/data/SunCross/SunRILs.yaml")

###Performing analysis
GWAS <- read.delim("output/run_2022-08-18/SunRILs_gwas_2022.csv", sep=",")
GWAS_subset <- data.frame(id = GWAS$id, p = GWAS$p_Awns)
GWAS_subset$chr <- str_replace(str_replace(GWAS_subset$id, "^S", ""), "_\\d*$", "")
GWAS_subset$pos <- as.numeric(str_replace(GWAS_subset$id, "^S\\d[ABD]_", ""))

chrom_colors = c("#5D8712", "#6B4300", "#2C397A")
SunValues <- data.frame()
families <- unique(genotype@ped$family)
for (j in families) {
	SunEffect <- data.frame()
	SunLOD <- data.frame()
	SunCross_sub <- SunCross[ind = ind_ids(SunCross)[grep(j, ind_ids(SunCross))]]
	#SunCross_sub <- SunCross
	dir <- paste("output/effects/", Sys.Date(), sep="")
	dir.create(dir)
	SunProb <- calc_genoprob(cross=SunCross_sub, map=NULL)
	SunKinship <- calc_kinship(SunProb)
	for (i in chr_names(SunCross_sub)) {
		##Do analysis
		##Get LOD score
		LOD <- scan1(SunProb[,i],SunCross_sub$pheno[,"Awns"], SunKinship)
		##Get effect estimate coefficients with BLUP scores
		SunEff <- scan1coef(SunProb[,i], SunCross_sub$pheno[,"Awns"], SunKinship)
		#SunBlup <- scan1blup(SunProb[,i], SunCross$pheno[,"Raleigh"], SunKinship)
		##Add to unified dataframe file
		SunEffect <- rbind(SunEffect, SunEff)
		SunLOD <- rbind(SunLOD, LOD)
	}
	##Save data in unified dataframe and format
	SunVal <- as.data.frame(merge(SunEffect, SunLOD, by=0))
	names(SunVal)[names(SunVal) == "pheno1"] <- "LOD" 
	names(SunVal)[names(SunVal) == "Row.names"] <- "id"
	SunVal$Family <- j
	SunValues <- rbind(SunValues, SunVal)
	SunVal <- merge(SunVal, data.frame(id = genotype@snps$id, chr = genotype@snps$chr, pos = genotype@snps$pos), by='id')
	###Plot out with GWAS for comparisoZ
	png(width=2500, height=1500, pointsize = 15, filename = paste(dir, "/", j, "_", 'effect.png', sep=""))
	#png(width=2500, height=1500, pointsize = 15, filename = paste(dir, "/", 'all_pop_awn_test.png', sep=""))
	par(mfrow=c(3,1), mar=numeric(4), oma = c(4, 4, .5, .5), mgp = c(2, .6, 0))
	layout(mat = matrix(c(1,2,3), nrow=3, ncol=1), heights = c(2,1,1))
	# par(mfrow=c(2,1), mar=numeric(4), oma = c(4, 4, .5, .5), mgp = c(2, .6, 0))
	# layout(mat = matrix(c(1,2), nrow=2, ncol=1), heights = c(2,1))

	#plot_coefCC(as.matrix(SunEffect), SunCross_sub$gmap, scan1_output=as.matrix(SunLOD), columns=1:2, main="Height Effect Scores", col=c("blue","red"))
	#plot_coef(as.matrix(SunEffect), SunCross$gmap, columns=1:2, main="Height Effect Scores", col=c("blue","pink"), xaxt="n")
	effect_df <-  data.frame(chr = SunVal$chr, pos = SunVal$pos, p = 10^SunVal$AA)
	manhattan(effect_df, cex=3, xaxt="n", chrom.col = chrom_colors)
	axis(2L)
	box()
	#axis(1L)
	#plot(SunLOD$pheno1, type='l', col = "purple", axes = FALSE, xaxt="n")#, SunCross$gmap[])
	LOD_df <- data.frame(chr = SunVal$chr, pos = SunVal$pos, p = 10^(-SunVal$LOD))
	manhattan(LOD_df, cex=1, axes=FALSE, xaxt="n", chrom.col = chrom_colors)
	#axis(1L)
	axis(2L)
	box()
	manhattan(GWAS_subset, cex=1, chrom.col = chrom_colors, axes=FALSE)#, xlim=c(min(GWAS_subset$id),max(GWAS_subset$id)))
	#axis(1L)
	axis(2L)
	box()
	dev.off()
}
write.csv(SunValues, paste(dir, "/", "SunRILs_effect_LOD.csv", sep=""))


SunProb <- calc_genoprob(cross=SunCross, map=NULL)
g <- maxmarg(SunProb, SunCross$gmap, chr='5A', pos=698528417, return_char=TRUE)
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, SunCross$pheno[,"Awns"], ylab="Awn phenotype")


Columns <- c("Marker", "Linkage Group (Chr)", "Location", "Trait", "Environment (Loc)", "LOD", "R2", "Additive Effect")
summary_table <- subset(SunVal, LOD >= sort(SunVal$LOD, decreasing=TRUE)[length(SunVal$LOD)*.001], select=c(id, chr, pos, LOD))




##Replot without largest effects in order to be able to see real resolution
reduced_SunEffect <- subset(SunEffect, abs(AA) < 50)
for (i in chr_names(SunCross)) {
	SunEff <- reduced_SunEffect[grep(i, rownames(reduced_SunEffect)),]
	comb <- merge(SunEff, SunLOD, by=0)
	rownames(comb) <- comb$Row.names
	comb <- as.matrix(comb[-1])
	png(width=2500, height=1500, pointsize = 15, filename = paste("output/effects/", i, '_effect.png', sep=""))
	plot_coefCC(comb[,1:3], SunCross$gmap[i], scan1_output=as.matrix(comb[,4]), columns=1:2, main=paste(i, "Effect Score", sep=""), col=c("blue", "pink"))
	dev.off()
}



g <- maxmarg(SunProb, SunCross$gmap, chr="6A", pos=6e+08, return_char=TRUE)
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, SunCross$pheno[,"Raleigh"], ylab="height phenotype 6A QTL spike")




par(mfrow=c(3,1), mar=numeric(4), oma = c(4, 4, .5, .5), mgp = c(2, .6, 0))
#plot_coefCC(as.matrix(SunEffect), SunCross$gmap, scan1_output=as.matrix(SunLOD), columns=1:2, main="Height Effect Scores", col=c("blue","pink"), ylim=c(0,100))
plot_coef(as.matrix(SunEffect), SunCross$gmap, columns=1:2, main="Height Effect Scores", col=c("blue","pink"), ylim=c(0,100), mar=c(0,1,1,1), axes = FALSE, xaxt="n")
#axis(2L)
box()
#plot(SunLOD, SunCross$gmap,mar=c(0,1,0,1))
manhattan(GWAS_subset, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), mar=c(0,1,0,1), axes = FALSE, xaxt="n")
axis(2L)
box()
manhattan(GWAS_subset, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), mar=c(1,1,0,1))
axis(1L)
axis(2L)
box()
