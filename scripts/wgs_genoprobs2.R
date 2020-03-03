#!/usr/bin/env Rscript

library("data.table")
library("tibble")
library("dplyr")

founders=c("A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra",
"DK63","EP1_inra","F492","FV252_inra","FV2_inra","ND245","OH43_inra","VA85","W117_inra")

#Take in chromosome 1..10 as arguments
args=commandArgs(trailingOnly=TRUE)
c = as.character(args[1])

#Read in WGS positions for chromosome
#wgs=fread(sprintf('biogemma/WGS_chr%s.txt',c),data.table=F,header=F)
#names(wgs)=c('chr','pos')
#wgs=wgs %>% mutate(marker=paste0('S',chr,'_',pos))
#wgs=wgs[,c('marker','chr','pos')]

options(scipen=999)

#for each chromosome
#read in files
pr=readRDS(sprintf('Biogemma_082318/600K_DHgenoprobs/bg%s_0823_genoprobs.rds',c))

#Read in physical map
#pmap=fread(sprintf('Biogemma_082318/Biogemma_pmap_c%s.csv',c),data.table=F)
#pmap$pos=pmap$pos*1e6

fgeno=fread(sprintf('biogemma/Biogemma_WGS_all_alleles_final_chr%s.txt',c),data.table=F)
#names(fgeno)=c('marker',founders)
fgeno %>% mutute(pos=strsplit(marker,'_')[[1]][2])
fgeno$pos=as.numeric(fgeno$pos)
fgeno=fgeno[,c('marker','pos',founders)]

#positions=rbind(fgeno[,c(,'marker','chr','pos')],pmap)
#positions=positions[order(positions$pos),]
#rownames(positions)=seq(1,nrow(positions))

#Read in founder WGS genotypes
#fgeno=fread(sprintf('Biogemma_080618/Biogemma_0904_WGS_foundergeno_chr%s.csv',c),data.table=F)
#rownames(fgeno)=fgeno[,1]
#reorder founders
#fgeno=fgeno[c('ind',founders),]



#len=dim(fgeno)[2]
#bin_fgeno=ifelse(fgeno[1:len,3:18]=='A',0,ifelse(fgeno[1:len,3:18]=='B',1,NA))
#Convert to 0 and 1

#remove(fgeno)

#Normalize probabilities based on dropped values

drop_missing <- function(alleles,probs){
    x=alleles[!is.na(alleles)]
    y=probs[!is.na(alleles)]*(1/sum(probs[!is.na(alleles)]))
    return(sum(x*y))
}


#Combine 600K and WGS positions together
#positions=rbind(wgs,pmap)
#positions=positions[order(positions$pos),]
#rownames(positions)=seq(1,nrow(positions))

#remove(wgs)
all_probs=c()
for( i in seq(1,344)){
    ind=pr[[1]][i,,]
    ind=t(ind)
    ind=as.data.frame(ind)
    names(ind)=c("marker",founders)
    ind=merge(ind,pmap,by.x="marker",by.y="marker")
    ind=ind[order(ind$pos)]
    rownames(ind)=seq(1,nrow(ind))
    ind=ind[,c("marker","chr","pos",founders)]

    f_probs=c()
    for(f in seq(1,16)){
        founder=founders[f]
	f_ind=ind[,c('pos',founder)]
	f_interp=approxfun(f_ind$pos,f_ind[,c(founder)],method="linear",yleft=f_ind[,c(founder)][1],yright=f_ind[,c(founder)][len])
	f_probs=rbind(wgs_f,f_interp(fgeno$pos))
    }
    allele_probs=sapply(seq(1,len),function(x) drop_missing(bin_fgeno[1:16,x],f_probs[1:16,x]))
    all_probs=rbind(all_probs,allele_probs)
}

fwrite(all_probs,sprintf('Biogemma_WGS_chr%s_allele_probs.txt',c),sep='\t',row.names=F,quote=F)