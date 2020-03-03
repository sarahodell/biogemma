#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')
#founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#ped=fread('Biogemma_Founders_600K.ped',data.table=F)
#
#ped=fread(sprintf('Biogemma_Founders_600K_chr%s.ped',c),data.table=F)
#ped$V1=founders
#ped$V2=founders
#rownames(ped)=ped$V1
#ped=ped[ped$V1!="MBS847",]
#v1=ped$V1
#ped$V1=paste0(ped$V1,"_",ped$V2)
#ped$V2=paste0(v1,"_",ped$V2)
#ped=ped[,-c(2)]
#fwrite(ped,sprintf('Biogemma_Founders_600K_chr%s.ped',c),row.names=F,col.names=F,quote=F,sep='\t')

gmap=fread(sprintf('Biogemma_Founders_600K_chr%s_gmap.map',c),data.table=F)
names(gmap)=c("marker","chr","cM")
map=fread(sprintf('Biogemma_Founders_600K_chr%s.map',c),data.table=F)
names(map)=c("chr","marker","cM","pos")
m=match(map$marker,gmap$marker)
map$cM=gmap$cM[m]
fwrite(map[,c("chr","marker","cM","pos")],sprintf('Biogemma_Founders_600K_chr%s.map',c),row.names=F,quote=F,col.names=F,sep='\t')

