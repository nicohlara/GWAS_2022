#####Using the qtl2 package, create a cross2 object and perform analysis on it
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-13

###SET WORKING DIRECTORY
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
library(tidyverse)
library(gaston)
library(qtl2)

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")

###creating cross2 files
#filter to significant markers
#sig_markers <- read.csv("output/SunRILs_sig_mark_2022.csv", sep=",")
#filter by height
#sig_markers <- subset(sig_markers, trait=="Height")# & p < 1e-10)
#geno_subset <- select.snps(genotype, id %in% sig_markers$marker)
#geno_subset <- select.inds(geno_subset, family %in% c("UX2029", "UX2000", "UX1995"))
#geno <- as.matrix(geno_subset)
geno <- as.matrix(genotype)
geno[geno==0] <- "AA"
geno[geno==1] <- "AB"
geno[geno==2] <- "BB"
write.csv(geno, "output/data/SunCross/SunRILs_geno.csv")
#gmap <- subset(geno_subset@snps, select=c("id", "chr", "pos"))
gmap <- subset(genotype@snps, select=c("id", "chr", "pos"))
names(gmap)[names(gmap)=="id"] <- "marker"
write.csv(gmap, "output/data/SunCross/SunRILs_gmap.csv", row.names=FALSE)
pheno <- subset(phenotype, select=c("Entry", "Location", "Height"))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995"))
pheno_K <- subset(pheno, Location == "Kinston", select=c("Entry", "Height"))
pheno_K <- rename(pheno_K, Kinson = Height)
pheno_R <- subset(pheno, Location == "Raleigh", select=c("Entry", "Height"))
pheno_R <- rename(pheno_R, Raleigh = Height)
pheno <- merge(pheno_R, pheno_K, by="Entry")
write.csv(pheno, "output/data/SunCross/SunRILs_pheno.csv", row.names=FALSE)
#phenocovar <- data.frame(pheno=c("Kinston", "Raleigh"), location = c("Kinston", "Raleigh"))
#write.csv(phenocovar, "/Users/nico/Desktop/cross2_file_test/SunRILs_phenocovar.csv", row.names=FALSE)
covar <- distinct(subset(phenotype, select=c("Entry", "Cross_ID")))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995")))
write.csv(covar, "output/data/SunCross/SunRILs_covar.csv", row.names=FALSE)
SunCross <- read_cross2("output/data/SunCross/SunRILs.yaml")

# ###Performing analysis
# pr <- calc_genoprob(cross=SunCross, map=NULL)
# c2eff_me <- scan1coef(pr[,"6A"], SunCross$pheno[,"Raleigh"])#iron$pheno[,"liver"])
# plot_coef(c2eff_me, map=SunCross$gmap, columns=1:2)
# #c2effB <- scan1coef(pr[,"6A"], SunCross$pheno[,"Raleigh"], contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
# # include BLUPs
# c2blup <- scan1blup(pr[,"6A"], SunCross$pheno[,"Raleigh"])
# plot_coef(c2eff_me, SunCross$gmap["6A"], columns=1:3)
# plot(c2blup, SunCross$gmap["6A"], columns=1:3, add=TRUE, lty=2, legend = "topright")

chrom_colors = c("#5D8712", "#6B4300", "#2C397A")
#SunEffect <- data.frame()
#SunLOD <- data.frame()
#SunCross_sub <- SunCross
SunValues <- data.frame()
families <- unique(genotype@ped$family)
for (j in families) {
	SunEffect <- data.frame()
	SunLOD <- data.frame()
	SunCross_sub <- SunCross[ind = ind_ids(SunCross)[grep(j, ind_ids(SunCross))]]
	dir <- paste("output/effects/", Sys.Date(), sep="")
	dir.create(dir)
	# dir.create("output/effects")
	SunProb <- calc_genoprob(cross=SunCross_sub, map=NULL)
	SunKinship <- calc_kinship(SunProb)
	for (i in chr_names(SunCross_sub)) {
		##Do analysis
		##Get LOD score
		LOD <- scan1(SunProb[,i],SunCross_sub$pheno[,"Raleigh"], SunKinship)
		##Get effect estimate coefficients with BLUP scores
		SunEff <- scan1coef(SunProb[,i], SunCross_sub$pheno[,"Raleigh"], SunKinship)
		#SunBlup <- scan1blup(SunProb[,i], SunCross$pheno[,"Raleigh"], SunKinship)
		##plot out and save
		#png(width=2500, height=1500, pointsize = 15, filename = paste(dir, "/", j, "_", i, "_", 'effect.png', sep=""))
		#plot_coefCC(SunEff, SunCross_sub$gmap[i], scan1_output=LOD, columns=1:2, col=c("blue", "pink"), main=paste(j, i, "Effect Score", sep=" "))
		#plot(SunBlup, SunCross$gmap[i], columns=1:2, add=TRUE, lty=2, legend = "topright")
		#dev.off()
		##Add to main files
		#cbind(SunEff, family = j)
		#SunLOD$family <- j
		SunEffect <- rbind(SunEffect, SunEff)
		SunLOD <- rbind(SunLOD, LOD)
	}
	
	#sort(SunEffect$AA)
	#largest_SunEffect <- subset(SunEffect, abs(AA) > 50)
	#write.csv(largest_SunEffect, "output/effects/SunRILs_height_largest_effects.csv")
	#write.csv(SunEffect, paste(dir, "/", j, "_", "SunRILs_height_effects.csv", sep=""))
	#write.csv(SunLOD, paste(dir, "/", j, "_", "SunRILs_LOD.csv", sep=""))
	SunVal <- as.data.frame(merge(SunEffect, SunLOD, by=0))
	names(SunVal)[names(SunVal) == "pheno1"] <- "LOD" 
	names(SunVal)[names(SunVal) == "Row.names"] <- "id"
	SunVal$Family <- j
	SunValues <- rbind(SunValues, SunVal)
	SunVal <- merge(SunVal, data.frame(id = genotype@snps$id, chr = genotype@snps$chr, pos = genotype@snps$pos), by='id')
 	
	
	###Plot out with GWAS
	png(width=2500, height=1500, pointsize = 15, filename = paste(dir, "/", j, "_", 'effect.png', sep=""))
	

	layout(mat = matrix(c(1,2,3), nrow=3, ncol=1), heights = c(2,1,1))
	 
	#par(mfrow=c(3,1), mar=numeric(4), oma = c(4, 4, .5, .5), mgp = c(2, .6, 0))
	#plot_coefCC(as.matrix(SunEffect), SunCross_sub$gmap, scan1_output=as.matrix(SunLOD), columns=1:2, main="Height Effect Scores", col=c("blue","red"))
	#plot_coef(as.matrix(SunEffect), SunCross$gmap, columns=1:2, main="Height Effect Scores", col=c("blue","pink"), xaxt="n")
	effect_df <-  data.frame(chr = SunVal$chr, pos = SunVal$pos, p = 10^SunVal$AA)
	manhattan(effect_df, cex=3, xaxt="n", chrom.col = chrom_colors)
	axis(2L)
	box()
	#axis(1L)
	#plot(SunLOD$pheno1, type='l', col = "purple", axes = FALSE, xaxt="n")#, SunCross$gmap[])
	LOD_df <- data.frame(chr = SunVal$chr, pos = SunVal$pos, p = SunVal$LOD)
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



GWAS <- read.delim("output/run_2022-08-18/SunRILs_gwas_2022.csv", sep=",")
GWAS_subset <- data.frame(id = GWAS$id, p = GWAS$p_Height)
GWAS_subset <- dataframe[complete.cases(dataframe),]
GWAS_subset <- GWAS_subset[complete.cases(dataframe),]
GWAS_subset <- GWAS_subset[complete.cases(GWAS_subset),]
GWAS_subset $LOG <- -log10(GWAS_subset$p)
GWAS_subset $LOG <- -log10(GWAS_subset$p)
GWAS_subset $chr <- str_replace(str_replace(GWAS_subset $id, "^S", ""), "_\\d*$", "")
GWAS_subset $pos <- as.numeric(str_replace(GWAS_subset $id, "^S\\d[ABD]_", ""))
GWAS_subset $FDR <- p.adjust(GWAS_subset $p, method = "fdr")
GWAS_subset $bon <- p.adjust(GWAS_subset $p, method = "bonferroni")

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
