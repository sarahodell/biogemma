a#!/usr/bin/env Rscript

library('data.table')
library('reshape2')
library('ggplot2')
library('tidyverse')

bounds=fread('../../../GridLMM/Biogemma_QTL.csv',data.table=T)
ftdays=c("male_flowering_days","female_flowering_days")
bounds=bounds[!(bounds$Phenotype %in% ftdays),]
rownames(bounds)=seq(1,nrow(bounds))
ld=fread('ld_bundled_links3.txt',data.table=F)
ld$V8=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(strsplit(ld$V7[x],',')[[1]][1],'=')[[1]][2]))
ld$V9=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(strsplit(ld$V7[x],',')[[1]][2],'=')[[1]][2]))
ld=ld[,-7]
names(ld)=c('chr1','start1','end1','chr2','start2','end2','size1','size2')
ld$chr1=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld$chr1[x],'r')[[1]][2]))
ld$chr2=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld$chr2[x],'r')[[1]][2]))

png('ld_region_sizes.png',width=800,height=800)
ggplot(ld,aes(x=size1/1e6,y=size2/1e6)) + geom_point() + xlab("Region Size (Mb)") + ylab("Region Size (Mb)") + ggtitle("Interchromosomal High LD (R-squared >= 0.8) Region Size")
dev.off()

bounds$pheno_env_id=paste0(bounds$pheno_env,'_',bounds$ID)
qtl_regions=data.table(bounds[,c('pheno_env_id','Phenotype','Environment','Method','ID','Chromosome','left_bound_bp','alt_right_bound_bp')])
#ld1=data.table(ld[,c('chr1','start1','end1')])
#ld2=data.table(ld[,c('chr2','start2','end2')])
ld=data.table(ld)
setkey(ld,chr1,start1,end1)
comparison1=foverlaps(qtl_regions,ld,by.x=c('Chromosome','left_bound_bp','alt_right_bound_bp'),by.y=c('chr1','start1','end1'),nomatch=NA)

setkey(ld,chr2,start2,end2)
comparison2=foverlaps(qtl_regions,ld,by.x=c('Chromosome','left_bound_bp','alt_right_bound_bp'),by.y=c('chr2','start2','end2'),nomatch=NA)
