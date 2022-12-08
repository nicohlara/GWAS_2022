#####Using the qtl2 package, perform CIM on a cross2 object
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-13

###SET WORKING DIRECTORY
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
library(tidyverse)
library(gaston)
library(qtl2)

##load in object
SunCross <- read_cross2("output/data/SunCross/SunRILs.yaml")
##insert pseudomarkers into genetic map and increase marker density
##probably unneccessary 
#map <- insert_pseudomarkers(SunCross$gmap, step=1)
map <- SunCross$gmap
##calculate QTL probabilities
pr <- calc_genoprob(SunCross, error_prob=0.002, cores=4)
##calculate a kinship matrix of relationship among individuals
kinship <- calc_kinship(pr)
##Perform genome scan using Haley-Knott regression
out <- scan1(pr, SunCross$pheno, kinship, cores=8)
##plot it
# ymx <- maxlod(out[,'WeightOf1000Particles'])
# plot(out, SunCross$gmap, lodcolumn ='WeightOf1000Particles', col='slateblue', ylim=c(0, ymx*1.2))
# plot(out, SunCross$gmap, lodcolumn ='SampleAreaAverage', col='violetred', add=TRUE)
# legend("topleft", lwd=2, col=c("slateblue", 'violetred'), c("Seed Weight", "Seed Size"), bg='gray90')
png(width = 2000, height = 1500, filename='CIM/output/raleigh_height_LOD.png')
plot(out, map, lodcolumn = 'Height')
dev.off()
png(width = 2000, height = 1500, filename='CIM/output/raleigh_ave_SpS_LOD.png')
plot(out, map, lodcolumn = 'ave_SpS')
dev.off()


##binary trait scan for awns



##obtain estimated coefficients
fivea_height <- scan1coef(pr[,"5A"], SunCross$pheno[,'Height'], kinship)
plot_coef(fivea_height, map, legend='topright')
summary(SunCross$pheno[,'Height'])

a <- find_peaks(out, map, threshold=4, drop=1.5)
g <- maxmarg(pr, map, chr='7B', pos=403036530, return_char=TRUE)
png(width = 1500, height = 1500, filename='CIM/output/raleigh_height_qtl_effect.png')
plot_pxg(g, SunCross$pheno[,'Height'], main='Largest QTL Effect for Height', sub="7B-403036530")
dev.off()
g <- maxmarg(pr, map, chr='1B', pos=6474600, return_char=TRUE)
png(width = 1500, height = 1500, filename='CIM/output/raleigh_sps_qtl_effect.png')
plot_pxg(g, SunCross$pheno[,'ave_SpS'], main='Largest QTL Effect for SpS', sub="1B-6474600")
dev.off()

# continue = FALSE
# for (i in rownames(a)) {
#   print(a[i, 'lodcolumn'])
#   g <- maxmarg(pr, map, chr=a[i,'chr'], pos=a[i,'pos'], return_char=TRUE)
#   plot_pxg(g, SunCross$pheno[,a[i, 'lodcolumn']])
#   cont_q <- readline('Enter any key to proceed')
# }
