###Analysis of 2022 Raleigh location height trait using CIM, testing methodology
##Created by Nicolas Lara
##Last modified 2022-12-5 

##in r
file_path <- "/project/guedira_seq_map/nico/r_packages"
.libPaths(file_path)
##Had to do manually, not sure why
#install.packages('tidyverse', repos='https://cran.case.edu/')
#install.packages('qtl', repos='https://cran.case.edu/')
library(tidyverse)
library(qtl)
setwd('/home/nicolas.lara/GWAS')


CIM_matrix <- read.cross(format = 'csv', file = 'CIM_matrix.csv', genotypes = c('1','2','0'))
scan.cim <- cim(CIM_matrix, pheno.col=2, map.function='kosambi')
scan.cim.perm = cim(CIM_matrix, pheno.col=2, map.function="kosambi", n.perm=1000)
write.csv(summary(scan.cim.perm), 'output/scan_cim_perm.csv')
write.csv(summary(scan.cim, threshold = 0.05), 'output/scan_cim.csv')
png(width=2500, height=1500, filename='output/scan_cim.png')
plot(scan.cim)
dev.off()
#png(width=2500, height=1500, filename='output/marker_effect.png')
#plotPXG(CIM_matrix, pheno.col=2)
#dev.off()


