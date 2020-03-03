#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

### Analyzing DH line genotype probabilities


library("qtl2")
library("dplyr")
library("data.table")
library("ggplot2")
library("reshape2")
library("tibble")


all_pr=readRDS(sprintf('bg%s_genoprobs_010319.rds',c))
all_pr_c=clean_genoprob(all_pr)
pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

hex_colors=c("#f42896","#84ef7c","#234489","#8ed1d6","#702349","#f2875b","#f1ff00","#bbbbbb","#ff2a3a","#56cc59","#663dd3","#478959","#47145b",'#0f0e0e',"#ad147a", "#afb735","#ff5a00","#fc1919")
founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra",
"FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

coverage=c()

for(f in 1:16){
    df=all_pr_c[[1]][,f,]
    tm =  t(df)
    mdf=as.data.frame(tm,row.names=rownames(tm))
    mdf=rownames_to_column(mdf,"marker")
    mdf = merge(mdf,pmap,by.x="marker",by.y="marker")
    mdf = mdf[,c(2:345,347)]
    mlong=melt(mdf,id="pos")
    max_prob=data.frame(pos=mdf$pos,total=rowSums(mdf[1:344]),founder=rep(founders[f],dim(mdf)[1]),chr=rep(c,dim(mdf)[1]))
    coverage=rbind(coverage,max_prob)
}
cov_counts=coverage %>% group_by(pos,founder) %>% summarize(tot_cov=round(total,2))

fwrite(cov_counts,sprintf('founder_coverage_counts_chr%s.txt',c),sep='\t',row.names=F,quote=F)

samples=names(all_pr[[1]][,1,1])
options(warn=-1)

xo_no=c()
all_breaks=c()
for (m in 1:344){
  md <- all_pr_c[[1]][m,,] #pull out the probs for one individual
  tm <- t(md) #transpose matrix
  mdf <- as.data.frame(tm,row.names = rownames(tm))
  mdf<-rownames_to_column(mdf,"marker")
  names(mdf)=c("marker",founders)
  mdf <- merge(mdf,pmap,by.x='marker',by.y='marker') #get the physical marker positions
  mdf <- mdf[,c(2:17,19)]
  mlong<-melt(mdf,id="pos") #melt the dataframe
  max_prob<-mlong %>% group_by(pos) %>% top_n(1,value) #grab the highest probability parent for each position
  max_prob<-as.data.frame(max_prob)
  max_prob=max_prob[max_prob$value >= 0.90,]
  max_prob<-max_prob[order(max_prob$pos),] #sort in ascending order
  rownames(max_prob)=seq(1:nrow(max_prob))
  # Find the most probable breakpoint between donors
  breaks=c()
  index=which(max_prob$variable != dplyr::lag(max_prob$variable))
  xo_no=c(xo_no,length(index))
  for(i in 1:length(index)){
    sample=samples[m]
    chr=10
    if(i==1){
      start=min(pmap$pos)
      end=max_prob[index[i]-1,'pos']
      donor=as.character(max_prob[1,'variable'])
      breaks=rbind(breaks,c(sample,chr,start,end,donor,donor))
    }
    start=max_prob[index[i],'pos']
    donor=as.character(max_prob[index[i],'variable'])
    prev=i-1
    if(i==length(index)){
      end=max(max_prob$pos)
    }
    else{
      end=max_prob[c(index[i+1]-1),'pos']
    }

    breaks=rbind(breaks,c(sample,chr,start,end,donor,donor))
  }
  breaks<-as.data.frame(breaks)
  names(breaks)=c('sample','chr','start','end','donor1','donor2')
  all_breaks=rbind(all_breaks,breaks)
}
options(warn=0)

fwrite(all_breaks,sprintf('Biogemma_breakpoints_chr%s.txt',c),row.names=F,quote=F,sep='\t')

xo=data.frame(sample=samples,chr=10,xo_no=xo_no)
fwrite(xo,sprintf('Biogemma_xo_number_chr%s.txt',c),row.names=F,quote=F,sep='\t')
