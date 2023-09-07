###Create some summary plots of phenotypic data
##Created by Nicolas A. H. Lara
##Last edit: 2022-12-6

library(tidyverse)
library(yarrr)

### SET WORKING DIRECTORY ###
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")

column_rename <- c('Awns', 'Winter Dormancy Release', 'Powdery Mildew', 'Septoria Nodorum Blotch', 'Height', 'Leaf Rust', 'Spikelets per Spike', 'Infertile Spikelets per Spike', 'Seed Weight', 'Seed Size (mean)', 'Seed Size (median)', 'Days to Head', 'Seeds per Spikelet')
###Correlation matrix of all variables
cor_matrix <- phenotype %>% 
	select(-c(Entry, Location, Cross_ID)) %>%
	'colnames<-'(column_rename) %>%
	as.matrix %>%
	cor(use = 'pairwise.complete.obs') %>%
	round(2) %>%
	as.data.frame %>%
	rownames_to_column(var = 'var1') #%>%
	gather(var2, value, -var1) %>%
	arrange(desc(value))

ggplot(data = cor_matrix, aes(x=var1, y=var2, fill=value)) + 
	geom_tile() + 
	scale_fill_gradient2(low='#CC3300', high='#009900', limit=c(-1,1), name="Correlation") +
	geom_text(aes(var2, var1, label = value), color = "white", size=3) +
	theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank())
ggsave('output/summary_plots/cor_heatmap.png')

pop_sum <- phenotype %>% select(c(Location, Cross_ID, Entry)) %>% unique
ggplot(pop_sum, aes(x=Cross_ID, fill = Location)) + geom_bar(position='dodge') + 
	theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())
ggsave('output/summary_plots/pop_counts.png')

png(width=2500, height=1500, filename='output/summary_plots/phenotype_pirate.png')
par(las=2, mfrow=c(1,2), cex=2)
pirateplot(Height ~ Cross_ID, data = subset(phenotype, Cross_ID != 'Parent'))
pirateplot(ave_SpS ~ Cross_ID, data = subset(phenotype, Cross_ID != 'Parent'))
dev.off()