#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')
library('tidyverse')

int_1=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_600K_germline_IBD_chr%s.bed',c),data.table=F,header=F)
names(int_1)=c('chr','start','end','name')
int_2=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_germline_IBD_chr%s.bed',c),data.table=F,header=F)
names(int_2)=c('chr','start','end','name')

for(n in unique(int_1$name)){
  fwrite(int_1[int_1$name==n,],sprintf('ibd_segments/comparison/tmp/%s_600K_c%s.bed',n,c),row.names=F,quote=F,col.names=F,sep='\t')
}

for(n in unique(int_2$name)){
  fwrite(int_2[int_2$name==n,],sprintf('ibd_segments/comparison/tmp/%s_wgs_c%s.bed',n,c),row.names=F,quote=F,col.names=F,sep='\t')
}


#int_join=inner_join(int_1,int_2)
#fwrite(int_join,sprintf('ibd_segments/comparison/Biogemma_GERMLINE_intersect_total_chr%s.bed',c),row.names=F,quote=F,sep='\t')

#Write summary file

#bed_wgs=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_germline_IBD_chr%s.bed',c),data.table=F)
#wgs_total=sum(bed_wgs$V3-bed_wgs$V2)


ibdseq=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_IBDSeq_chr%s.bed',c),data.table=F,header=F)
names(ibdseq)=c('chr','start','end','name')

for(n in unique(ibdseq$name)){
  fwrite(ibdseq[ibdseq$name==n,],sprintf('ibd_segments/comparison/tmp/%s_wgs_ibdseq_c%s.bed',n,c),row.names=F,quote=F,col.names=F,sep='\t')
}

refined=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_RefinedIBD_chr%s.bed',c),data.table=F,header=F)
names(refined)=c('chr','start','end','name')

for(n in unique(refined$name)){
  fwrite(refined[refined$name==n,],sprintf('ibd_segments/comparison/tmp/%s_RefinedIBD_c%s.bed',n,c),row.names=F,quote=F,col.names=F,sep='\t')
}
#Write out pair lists

# GERMLINE 600K and WGS
pair_list1=data.frame(name=intersect(int_1$name,int_2$name),stringsAsFactors=F)
fwrite(pair_list1,sprintf("ibd_segments/comparison/c%s_germline_pair_list.txt",c),row.names=F,quote=F,col.names=F,sep='\t')

pair_list2=data.frame(name=intersect(int_1$name,ibdseq$name),stringsAsFactors=F)
fwrite(pair_list2,sprintf("ibd_segments/comparison/c%s_600K_ibdseq_pair_list.txt",c),row.names=F,quote=F,col.names=F,sep='\t')

pair_list3=data.frame(name=intersect(int_2$name,ibdseq$name),stringsAsFactors=F)
fwrite(pair_list3,sprintf("ibd_segments/comparison/c%s_germline_wgs_ibdseq_pair_list.txt",c),row.names=F,quote=F,col.names=F,sep='\t')

# GERMLINE 600K and RefinedIBD

pair_list4=data.frame(name=intersect(int_1$name,refined$name),stringsAsFactors=F)
fwrite(pair_list4,sprintf("ibd_segments/comparison/c%s_germline_WGS_RefinedIBD_pair_list.txt",c),row.names=F,quote=F,col.names=F,sep='\t')
