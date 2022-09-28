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

