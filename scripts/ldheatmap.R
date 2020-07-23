#!/usr/bin/env Rscript

library('LDheatmap')
library('data.table')
library('tidyverse')


ld=fread('stats/ld_decay/Biogemma_DHLines_rsquared_all_chroms.ld',data.table=F)
map=fread('genotypes/plink_files/600K/Biogemma_Founders_600K.map',data.table=F)

ld_matrix = ld[,c('SNP_A','SNP_B','R2')] %>%
  group_by(`SNP_A`) %>%
  mutate(id = row_number()) %>%
pivot_wider(names_from='SNP_A',values_from='R2')
ld_matrix=as.data.frame(ld_matrix,stringsAsFactors=F)
rownames(ld_matrix)=ld_matrix$SNP_B
map=map[map$V2 %in% rownames(ld_matrix),]
ld_matrix=ld_matrix[,-1]
ld_matrix=as.matrix(ld_matrix)


png('stats/ld_decay/Biogemma_Founders_LD_heatmap.png')
LDheatmap(ld_matrix,genetic.distances=map$V3,distances=map$V4)
dev.off()
