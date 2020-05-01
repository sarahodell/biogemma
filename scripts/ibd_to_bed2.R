#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')

refined=fread(sprintf('ibd_segments/refinedibd/Biogemma_WGS_Founders_RefinedIBD_chr%s.ibd',c),data.table=F,header=F)
names(refined)=c('ID_1','ID_1_index','ID_2','ID_2_index','Chromosome','start','end','LOD','cM')
refined$pair=paste0(refined$ID_1,'-',refined$ID_2)
bed=refined[,c('Chromosome','start','end','pair')]

fwrite(bed,sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_RefinedIBD_chr%s.bed',c),row.names=F,quote=F,sep='\t',col.names=F)
