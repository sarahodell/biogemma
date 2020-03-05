#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')
library('tidyverse')

int_1=fread(sprintf('ibd_segments/comparison/Biogemma_GERMLINE_WGS_600K_chr%s_intersect.bed',c),data.table=F)
names(int_1)=c('chr','start','end','name')
int_2=fread(sprintf('ibd_segments/comparison/Biogemma_GERMLINE_600K_WGS_chr%s_intersect.bed',c),data.table=F)
names(int_2)=c('chr','start','end','name')


int_join=inner_join(int_1,int_2)
fwrite(int_join,sprintf('ibd_segments/comparison/Biogemma_GERMLINE_intersect_total_chr%s.bed',c),row.names=F,quote=F,sep='\t')

#Write summary file

bed_wgs=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_germline_IBD_chr%s.bed',c),data.table=F)
wgs_total=sum(bed_wgs$V3-bed_wgs$V2)
