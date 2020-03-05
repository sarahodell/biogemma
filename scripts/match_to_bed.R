#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library("data.table")
library("readr")

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

founders2=c("A632","B73","CO255","FV252","OH43", "A654","FV2","C103","EP1","D105","W117","B96","DK63","F492","ND245","VA85")

ibd_600K=read_table2(sprintf('ibd_segments/germline/600K/Biogemma_Founders_germline_IBD_chr%s.match',c),col_names=F)
ibd_600K=as.data.frame(ibd_600K,stringsAsFactors=F)
for(i in seq(1,dim(ibd_600K)[1])){
      if(!(ibd_600K[i,]$X1 %in% founders)){
        f2=which(founders2==ibd_600K[i,]$X1)
        ibd_600K[i,]$X1=founders[f2]
  
        f2=which(founders2==ibd_600K[i,]$X3)
        ibd_600K[i,]$X3=founders[f2]
  }
}

names(ibd_600K)=c('ID_1','Family_ID_1','ID_2','Family_ID2','Chromosome','left_pos','right_pos','left_marker','right_marker','n_markers','Genetic_length','Genetic_length_units','n_mismatch','ID1_Is_Homozygous','ID2_Is_Homozygous')

bed_600K=ibd_600K[,c('Chromosome','left_pos','right_pos')]
bed_600K$name=paste0(ibd_600K$ID_1,'-',ibd_600K$ID_2)
#bed_600K$score='.'
#bed_600K$strand=paste0(ibd_600K$ID_1,'-',ibd_600K$ID_2)

fwrite(bed_600K,sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_600K_germline_IBD_chr%s.bed',c),row.names=F,col.names=F,quote=F,sep='\t')

ibd_wgs=read_table2(sprintf('ibd_segments/germline/WGS/Biogemma_Founders_WGS_germline_IBD_chr%s.match',c),col_names=F)
ibd_wgs=as.data.frame(ibd_wgs,stringsAsFactors=F)

names(ibd_wgs)=c('ID_1','Family_ID_1','ID_2','Family_ID2','Chromosome','left_pos','right_pos','left_mar\
ker','right_marker','n_markers','Genetic_length','Genetic_length_units','n_mismatch','ID1_Is_Homozygous','ID2_Is_Homozygous')

bed_wgs=ibd_wgs[,c('Chromosome','left_pos','right_pos')]
bed_wgs$name=paste0(ibd_wgs$ID_1,'-',ibd_wgs$ID_2)
#bed_wgs$score='.'
#bed_wgs$strand=paste0(ibd_wgs$ID_1,'-',ibd_wgs$ID_2)

fwrite(bed_wgs,sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_germline_IBD_chr%s.bed',c),row.names=F,col.names=F,quote=F,sep='\t')
