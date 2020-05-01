#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')

ibdseq=fread(sprintf('ibd_segments/IBDSeq/Biogemma_WGS_Founders_IBDSEQ_chr%s.ibd',c),data.table=F,header=F)
names(ibdseq)=c('ID_1','ID_1_index','ID_2','ID_2_index','Chromosome','start','end','LOD')
ibdseq$pair=paste0(ibdseq$ID_1,'-',ibdseq$ID_2)
bed=ibdseq[,c('Chromosome','start','end','pair')]

fwrite(bed,sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_IBDSeq_chr%s.bed',c),row.names=F,quote=F,sep='\t',col.names=F)
