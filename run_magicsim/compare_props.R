#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
rep=as.character(args[[1]])
c=as.numeric(args[[2]])
cores=as.numeric(args[[3]])

library('data.table')
library('reshape2')
#library('tidyverse')
library('parallel')
library('MASS')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra",
"OH43_inra","A654_inra","FV2_inra","C103_inra",
"EP1_inra","D105_inra","W117_inra","B96","DK63",
"F492","ND245","VA85")

breaks=fread(sprintf('breaktables/MAGIC_DH_Sim_rep%s_breaktable.txt',rep))

pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
fprobs=readRDS(sprintf('qtl2_files/MAGIC_DHSim_rep%s_c%.0f_genoprobs.rds',rep,c))
dimnames(fprobs[[1]])[[2]]=founders

nind=dim(fprobs[[1]])[1]
nmarkers=dim(fprobs[[1]])[3]

actual_predicted=function(ind){
  indprobs=fprobs[[1]][ind,,]
  mdf=apply(indprobs,MARGIN=2,function(i) names(which.max(i)))
  #mdf=as.data.frame(matrix(mdf, ncol=400, byrow=FALSE),stringsAsFactors=F)
  #rownames(mdf)=dimnames(fprobs[[1]])[[3]]
  #names(mdf)=ind
  mdf=data.frame(marker=names(mdf),predicted=mdf,stringsAsFactors=F)
  mdf = merge(mdf,pmap,by.x="marker",by.y="marker")
  mdf=mdf[order(mdf$pos),]
  rownames(mdf)=mdf$marker
  #return(mdf)

  indbreaks=breaks[breaks$sample==ind & breaks$chr==c,]
  mdf$actual=sapply(seq(1,nmarkers), function(i) indbreaks[between(mdf$pos[i],indbreaks$start,indbreaks$end),]$donor)
  #rownames(t)=dimnames(fprobs[[1]])[[3]]
  #colnames(t)=dimnames(fprobs[[1]])[[1]]
  return(sum(mdf$actual==mdf$predicted)/nrow(mdf) * 100)
}

n_inds=paste0('Sim',seq(1,400))

print(system.time({
results=mclapply(n_inds,actual_predicted,mc.cores=cores)
}))
perc_acc=mean(unlist(results))
line=data.table(rep=rep,chr=as.numeric(c),perc_acc=perc_acc)
fwrite(line,file="genoprob_accuracy.txt",sep='\t',append=T,col.names=F,quote=F,row.names=F)


#mdf=sapply(seq(1,nind),function(x) sapply(seq(1,nmarkers),function(i) names(which.max(fprobs[[1]][x,,i]))))
#mdf=as.data.frame(matrix(mdf, ncol=400, byrow=FALSE),stringsAsFactors=F)
#rownames(mdf)=dimnames(fprobs[[1]])[[3]]
#colnames(mdf)=dimnames(fprobs[[1]])[[1]]
#mdf=rownames_to_column(mdf,"marker")
#mdf = merge(mdf,pmap,by.x="marker",by.y="marker")
#mdf=mdf[order(mdf$pos),]
#rownames(mdf)=seq(1,nrow(mdf))


#t=sapply(seq(1,400), function(x) sapply(seq(1,nmarkers), function(i) breaks[breaks$sample==paste0('Sim',x) & breaks$start<= mdf$pos[i] & breaks$end>mdf$pos[i] & breaks$chr==c,]$donor))


#match=sapply(seq(1,400), function(x) sum(mdf[,paste0('Sim',x)] == t[,paste0('Sim',x)]))
#perc_acc=mean(match/nmarkers)*100
#line=data.table(rep=rep,chr=as.numeric(c),perc_acc=perc_acc)
#fwrite(line,file="genoprob_accuracy.txt",sep='\t',append=T,col.names=F,quote=F,row.names=F)


#break_probs=c()
#for(n in seq(1,nind)){
#  ind=dimnames(fprobs[[1]])[[1]][n]
#  md=fprobs[[1]][n,,]
#  tm=t(md)
#  rownames(tm)=dimnames(fprobs[[1]])[[3]]
#  mdf=as.data.frame(tm,row.names=rownames(tm))
#  mdf=rownames_to_column(mdf,"marker")
#  mdf = merge(mdf,pmap,by.x="marker",by.y="marker")
  #mdf = mdf[,c(2:17,19)]
#  fmax=founders[apply(mdf[,2:17],MARGIN=1,which.max)]
#  max_prob=data.frame(ind=ind,pos=mdf$pos,founder=fmax,chr=c,marker=mdf$marker)
#  max_prob=max_prob[order(max_prob$pos),]
#  rownames(max_prob)=seq(1,nrow(max_prob))
#  break_probs=rbind(break_probs,max_prob)
#}
