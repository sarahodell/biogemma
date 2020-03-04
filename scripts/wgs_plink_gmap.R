#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')
founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#ped=fread('Biogemma_Founders_600K.ped',data.table=F)

ped=fread(sprintf('/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr%s.ped',c),data.table=F)
ped$V1=founders
ped$V2=founders
#rownames(ped)=ped$V1
#ped=ped[ped$V1!="MBS847",]
#v1=ped$V1
#ped$V1=paste0(ped$V1,"_",ped$V2)
#ped$V2=paste0(v1,"_",ped$V2)
#ped=ped[,-c(2)]
fwrite(ped,sprintf('/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr%s.ped',c),row.names=F,quote=F,sep='\t',col.names=F)

ogutmap=fread('/home/sodell/projects/biogemma/genotypes/misc/ogutmap_v4_ordered.csv',data.table=F)
gmap=ogutmap[ogutmap$chr==as.integer(c),]
approx_cM<-approxfun(gmap$pos,gmap$scaled_cM,yleft=0,yright=max(gmap$scaled_cM))
map=fread(sprintf('/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr%s.map',c),data.table=F)
names(map)=c("chr","marker","cM","pos")
map$cM=approx_cM(map$pos)
fwrite(map[,c("chr","marker","cM","pos")],sprintf('/home/sodell/projects/biogemma/genotypes/plink_files/WGS/Biogemma_Founders_WGS_chr%s.map',c),row.names=F,quote=F,col.names=F,sep='\t')

