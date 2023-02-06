###Data analysis playground
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-13

library(tidyverse)
library(gaston)
#library(lmerTest)
#library(RAINBOWR)
#source("analysis/GWAS_functions.R")
#install.packages('asreml')
library(asreml)
library(AGHmatrix)
library(ASRgenomics)
#install.packages("qtl2")
library(qtl2)
#install.packages('qtl')
library(qtl)

### SET WORKING DIRECTORY ###
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")


phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")
phenotype <- filter(phenotype, Entry %in% genotype@ped$id)


#t_ph_group is replicated averaged
#t_ph is not

t_ph_s1 <- subset(t_ph, Cross_ID != "Parent")

for (j in unique(t_ph_s1$Location)) {
	print(j)
	phen <- subset(t_ph_s1, Location == j)
	print(length(phen$Entry))
	print(length(unique(phen$Entry)))
	for (i in unique(phen$Cross_ID)) {
		print(i)
		print(length((subset(phen, Cross_ID == i)$Entry)))
		print(length(unique(subset(phen, Cross_ID == i)$Entry)))
	}
}
# for (i in unique(t_ph$Entry)) {
	# a <- subset(t_ph, Entry == i, select=c("Location", "Cross_ID", "Entry"))
	# if (length(a$Entry) < 2) { print(a) }
# }



###################Vrn segregation amid selected populations##################
geno_sub <- select.inds(genotype, family %in% c("UX1993", "UX2023", "UX1444", "Parents"))
geno_sub <- select.snps(geno_sub, (chr == "5A" & pos >= 585018041 & pos <= 598072246) | (chr == "5B" & pos >= 564000000 & pos <= 568000000))
matrix_geno <- as.matrix(geno_sub)
write.csv(matrix_geno, "/Users/nico/Documents/GitHub/GWAS_2022/output/vrn_gene_selection.csv")
geno_sub <- select.snps(geno_sub, id %in% c("S5A_585018041", "S5A_598072246", "S5B_564112036"))
write.csv(as.matrix(geno_sub), "/Users/nico/Documents/GitHub/GWAS_2022/output/vrn_gene_finalists.csv")



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
raw_pheno <- read.delim("output/data/2022_raw_phenotype.csv", sep=",")
##make genotype data for model
#rgeno <- t(as.matrix(genotype))
#marker <- str_split(rownames(rgeno), "_")
#chr <- as.matrix(lapply(marker, '[[', 1)) %>% 'colnames<-'('chrom')
#loc <- as.matrix(as.numeric(lapply(marker, '[[', 2))) %>% 'colnames<-'('pos')
#rgeno <- cbind(chr, loc, rgeno)
#geno <- as.matrix(genotype)
##The G matrix we calculated above is not positive definite. Use nearPD function from Matrix
#package to get a positive definite matrix and remove singularity.
##This comment taken from For726, may not be accurate
#G <- Gmatrix(geno, method = "VanRaden", ploidy = 2)
#RealizedPD = nearPD(G, keepDiag = T)
#G = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
#G = G + diag(0.01, nrow(G))
#attr(G, "dimnames") = RealizedPD[[1]]@Dimnames
#summary(eigen(G)$values)
Ginv.sparse <- G.inverse(G = G, sparseform = TRUE)$Ginv
#write.csv(Ginv.sparse, 'output/data/Ginv.csv', row.names=FALSE)
Ginv.sparse <- as.matrix(read.delim('output/data/Ginv.csv', sep=","))
rownames(Ginv.sparse) <- c(1:dim(Ginv.sparse)[1])
##Coerce factors
phenotype$Location <- as.factor(phenotype$Location)
phenotype$Entry <- as.factor(phenotype$Entry)
#make model
model1 <- asreml(Height ~ Location,
                 random = ~vm(Entry, Ginv.sparse),
                 data=phenotype)

gblup_mod <- asreml(Plant_Height ~ Environment,
                    random = ~vm(Genotype, Ginv.sparse),
                    data=rpheno)

#plot
plot(model1)


###CIM code
Sun_cross <- subset(phenotype, Location=='Raleigh', select=c('Entry', 'Awns'))
rownames(Sun_cross) <- Sun_cross$Entry
Sun_cross <- select(Sun_cross, -'Entry')
geno_matrix <- data.frame(as.matrix(genotype))
genotype_matrix <- merge(Sun_cross, geno_matrix, by=0, all=FALSE)
row.names(genotype_matrix) <- genotype_matrix$Row.names
genotype_matrix <- genotype_matrix[,-1]
col_names <- str_split(colnames(geno_matrix), "_")
chrom <- c(NA, sapply(col_names,"[[",1))
pos <- c(NA, as.numeric(sapply(col_names,"[[",2))/1e6)
genotype_matrix <- rbind(chrom, pos, genotype_matrix)

write.csv(genotype_matrix, "output/data/CIM_matrix.csv")
