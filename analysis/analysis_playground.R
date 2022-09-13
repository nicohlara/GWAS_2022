###Data analysis playground
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-13

library(tidyverse)
library(gaston)
#library(lmerTest)
#library(RAINBOWR)
#source("analysis/GWAS_functions.R")
#install.packages('asreml')
#library(asreml)
#install.packages("qtl2")
library(qtl2)


### SET WORKING DIRECTORY ###
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")


phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")
phenotype <- filter(phenotype, Entry %in% genotype@ped$id)

GWAS_loop <- function(genotype, phenotype) {
	plot_df <- data.frame()
	len <- length(genotype@snps$id)
	count <- 0
	print(Sys.time())
	for (i in c(1:len)) {
		if (round((i/len)*100) > count) {
			count <- round((i/len)*100)
			print(paste(count, "%", sep=""))
		}	
		a <- as.matrix(genotype[,i])
		b <- data.frame(Entry = rownames(a), Marker = as.vector(a[,1]))
		c <- colnames(a)
		if (length(unique(b$Marker)) > 1) {
			d <- merge(phenotype,b, by="Entry")
			p_vals <- data.frame(id = c)
			for (j in c(4:length(phenotype[1,]))) {
				if (length(unique(d$Location)) > 1) {
					###TREATING LOCATION AND FAMILY RANDOM EFFECTS
					mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Location) + (1|Cross_ID)))
					mm <- suppressMessages(lmer(data=d, ))
				} else {
					###Running without location
					mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Cross_ID)))
				}
				p <- anova(mm)["Marker",6]
				p_vals[paste("p_", names(phenotype[j]), sep="")] <- p
			}
			plot_df <- rbind(plot_df, p_vals)
			#plot_df <- rbind(plot_df, data.frame(id = c, p = p_val))
		}
	}
	return(plot_df)
}


height_pheno <- subset(phenotype, select=c("Location", "Cross_ID", "Entry", "Height"))
plot_df <- GWAS_loop(genotype, height_pheno)
png(width=2500, height=1500, pointsize = 15, filename = "output/test_plot")
sig_m <- manhattan_plot("Height test", plot_df$id, plot_df[,2])
dev.off()



###Experimenting with qtl2 package
setwd("/Users/nico/Desktop/mapping")
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
options(timeout=900) # set the download timeout to 900 seconds from the default 60 seconds to help with large file downloads
# these four commands download files from a url and places them in the data directory you created
download.file(url="https://ndownloader.figshare.com/files/18533342", destfile="./data/cc_variants.sqlite") 
download.file(url="https://ndownloader.figshare.com/files/24607961", destfile="./data/mouse_genes.sqlite")
download.file(url="https://ndownloader.figshare.com/files/24607970", destfile="./data/mouse_genes_mgi.sqlite")
download.file(url="ftp://ftp.jax.org/dgatti/qtl2_workshop/qtl2_demo.Rdata", destfile="./data/qtl2_demo.Rdata")
options(timeout=60) # reset the download timeout to the default 60 seconds
# create a function to query the SNP file, then use this new function  
# to select SNPs on chromosome 1 from 10 to 11 Mbp
snp_func = create_variant_query_func(dbfile = "~/Desktop/mapping/data/cc_variants.sqlite") 
snps = snp_func(chr = 1, start = 10, end = 11) 
# check the dimensions of this sample of the SNP file
dim(snps) 
# create a function to query the gene file, then select genes in the same region as above
gene_func = create_gene_query_func(dbfile = "~/Desktop/mapping/data/mouse_genes_mgi.sqlite") 
genes = gene_func(chr = 1, start = 10, end = 11) 
dim(genes) # check the dimensions


##What to copy
#read in phenotype file
iron <- read_cross2(file = system.file("extdata", "iron.zip", package="qtl2"))
#create genetic map, insert pseudomarkers
map <- insert_pseudomarkers(map=iron$gmap, step=1)
#calculate genotype probabilities of markers
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
#calculate QTL effects
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])
#plot coefficients
plot_coef(c2eff, map, legend = "topright")
#include additive and dominance effects for QTL scores
c2effB <- scan1coef(pr[,"2"], iron$pheno[,"liver"],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))
#include colors for plot
col <- c("green", "red")
plot_coef(c2effB, map["2"], columns=2:3, col=col)



###replicate with my data
pheno <- split(phenotype$Entry, phenotype$Heigh)
geno <- select.snps(genotype, chr == "1A")
geno <- data.frame(id = genotype@snps$id, MM = 1-genotype@snps$maf, Mm = 0, mm = genotype@snps$maf)
geno <- split(geno[-1], geno[1])
c2eff_me <- scan1coef(geno, pheno)



###creating cross2 files
#filter to significant markers
sig_markers <- read.csv("/Users/nico/Documents/GitHub/GWAS_2022/output/SunRILs_sig_mark_2022.csv", sep=",")
#filter by height
sig_markers <- subset(sig_markers, trait=="Height")# & p < 1e-10)
geno_subset <- select.snps(genotype, id %in% sig_markers$marker)
#geno_subset <- select.inds(geno_subset, family %in% c("UX2029", "UX2000", "UX1995"))
geno <- as.matrix(geno_subset)
geno[geno==0] <- "AA"
geno[geno==1] <- "AB"
geno[geno==2] <- "BB"
write.csv(geno, "/Users/nico/Desktop/cross2_file_test/SunRILs_geno.csv")
gmap <- subset(geno_subset@snps, select=c("id", "chr", "pos"))
names(gmap)[names(gmap)=="id"] <- "marker"
write.csv(gmap, "/Users/nico/Desktop/cross2_file_test/SunRILs_gmap.csv", row.names=FALSE)
pheno <- subset(phenotype, select=c("Entry", "Location", "Height"))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995"))
pheno_K <- subset(pheno, Location == "Kinston", select=c("Entry", "Height"))
pheno_K <- rename(pheno_K, Kinson = Height)
pheno_R <- subset(pheno, Location == "Raleigh", select=c("Entry", "Height"))
pheno_R <- rename(pheno_R, Raleigh = Height)
pheno <- merge(pheno_R, pheno_K, by="Entry")
write.csv(pheno, "/Users/nico/Desktop/cross2_file_test/SunRILs_pheno.csv", row.names=FALSE)
#phenocovar <- data.frame(pheno=c("Kinston", "Raleigh"), location = c("Kinston", "Raleigh"))
#write.csv(phenocovar, "/Users/nico/Desktop/cross2_file_test/SunRILs_phenocovar.csv", row.names=FALSE)
covar <- distinct(subset(phenotype, select=c("Entry", "Cross_ID")))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995")))
write.csv(covar, "/Users/nico/Desktop/cross2_file_test/SunRILs_covar.csv", row.names=FALSE)
SunCross <- read_cross2("/Users/nico/Desktop/cross2_file_test/SunRILs.yaml")

###Performing analysis
#map_me <- insert_pseudomarkers(map=SunCross$gmap, step=1)
pr <- calc_genoprob(cross=SunCross, map=NULL)
c2eff_me <- scan1coef(pr[,"6A"], SunCross$pheno[,"Raleigh"])#iron$pheno[,"liver"])
plot_coef(c2eff_me, map=SunCross$gmap)
#c2effB <- scan1coef(pr[,"6A"], SunCross$pheno[,"Raleigh"], contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
# include BLUPs
c2blup <- scan1blup(pr[,"6A"], SunCross$pheno[,"Raleigh"])
plot_coef(c2eff_me, SunCross$gmap["6A"], columns=1:3)
plot(c2blup, SunCross$gmap["6A"], columns=1:3, add=TRUE, lty=2, legend = "topright")









See(phenotype)
See(genotype)
pheno_subset <- subset(phenotype, Location=="Raleigh")

GWAS_pheno <- subset(pheno_subset, select=c("Entry", "Height"))
#covariate_factors <- subset(phenotype, select=c("Location"))
GWAS_geno <- data.frame(marker=genotype@snps$id, chrom=genotype@snps$chr, pos=genotype@snps$pos)
#rownames(GWAS_geno) <- GWAS_geno$marker
for (i in genotype@ped$id) {
	a <- as.data.frame(t(as.matrix(select.inds(genotype, id==i))))
	a$marker <- rownames(a)
	if (length(a[,1])>10) {
		GWAS_geno <- merge(GWAS_geno, a, by='marker', all.x=FALSE)
	}
}
GWAS_geno$chrom <- gsub("A", "1", GWAS_geno$chrom)
GWAS_geno$chrom <- gsub("B", "2", GWAS_geno$chrom)
GWAS_geno$chrom <- gsub("D", "3", GWAS_geno$chrom)
GWAS_geno$chrom <- as.numeric(GWAS_geno$chrom)
multi_snp <- RGWAS.multisnp(GWAS_pheno, GWAS_geno, haplotype=TRUE)# covariate.factor=covariate_factors, haplotype=TRUE)





###Gaston test
x <- genotype
K <- GRM(x)
eiK <- eigen(K)

sub_pheno <- subset(phenotype, Entry %in% genotype@ped$id)
y <- as.matrix(subset(sub_pheno, Location=='Raleigh', select=Height))

TAU <- seq(.5, 1.5, length=30)
S2 <- seq(1,3,length=30)
lik1 <- lmm.diago.likelihood(tau = TAU, s2 = S2, Y=y, eigenK = eiK)
H2 <- seq(0,1,length=51)
lik2 <- lmm.diago.likelihood(h2=H2, Y=y, eigenK = eiK)


Y <- y
X = matrix(1, nrow = length(Y))
mixed_model <- lmm.aireml(y, X, K)




lmm.aireml(Y, X = matrix(1, nrow = length(Y)), K,
EMsteps = 0L, EMsteps_fail = 1L, EM_alpha = 1,
min_tau, min_s2 = 1e-06, theta, constraint = TRUE, max_iter = 50L,
eps = 1e-05, verbose = getOption("gaston.verbose", TRUE),
contrast = FALSE, get.P = FALSE)



# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
# Compute Genetic Relationship Matrix
K <- GRM(x)
# eigen decomposition of K
eiK <- eigen(K)
# simulate a phenotype
set.seed(1)
y <- 1 + lmm.simu(tau = 1, sigma2 = 2, eigenK = eiK)$y
# Likelihood
TAU <- seq(0.5,1.5,length=30)
S2 <- seq(1,3,length=30)
lik1 <- lmm.diago.likelihood(tau = TAU, s2 = S2, Y = y, eigenK = eiK)
H2 <- seq(0,1,length=51)
lik2 <- lmm.diago.likelihood(h2 = H2, Y = y, eigenK = eiK)
# Plotting
par(mfrow=c(1,2))
lik.contour(TAU, S2, lik1, heat = TRUE, xlab = "tau", ylab = "sigma^2")
lines(lik2$tau, lik2$sigma2)
plot(H2, exp(lik2$likelihood), type="l", xlab="h^2", ylab = "likelihood")



###asreml
#get inverse of relationship matrix
gryphonped <- read.csv("/Users/nico/Documents/GitHub/wam_tuto/data/gryphonped.csv")
gryphon <- read.csv("/Users/nico/Documents/GitHub/wam_tuto/data/gryphon.csv")
ainv <- ainverse(gryphonped)
#make model
model1 <- asreml(
  fixed = bwt ~ 1, random = ~ vm(animal, ainv),
  residual = ~ idv(units),
  data = gryphon,
  na.action = na.method(x = "omit", y = "omit")
)
#plot
plot(model1)

