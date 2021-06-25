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
library('cowplot')

all_pr=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',c))
#all_pr=readRDS(sprintf('bg%s_genoprobs_010319.rds',c))
#all_pr_c=clean_genoprob(all_pr)
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

#hex_colors=c("#f42896","#84ef7c","#234489","#8ed1d6","#702349","#f2875b","#f1ff00","#bbbbbb","#ff2a3a","#56cc59","#663dd3","#478959","#47145b",'#0f0e0e',"#ad147a", "#afb735","#ff5a00","#fc1919")
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra",
"FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#coverage=c()

#all_highp=0
#for(c in 1:10){
#  all_pr=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%.0f_genoprobs.rds',c))
#  a=dim(all_pr[[1]])[1]
#  b=dim(all_pr[[1]])[3]
#  total=a*b
#  highp=0
#  per_ind=c()
#  thresh=0.8
#  for(f in 1:a){
#      df=all_pr[[1]][f,,]
#      r=sum(apply(df,MARGIN=2,function(x) length(x[x>=thresh])))
#      highp=highp+r
#      per_ind=c(per_ind,r/ncol(df))
#  }
#  all_highp=c(all_highp,highp/total)
#}



#for(f in 1:16){
#    df=all_pr_c[[1]][,f,]
#    tm =  t(df)
#    mdf=as.data.frame(tm,row.names=rownames(tm))
#    mdf=rownames_to_column(mdf,"marker")
#    mdf = merge(mdf,pmap,by.x="marker",by.y="marker")
#    mdf = mdf[,c(2:345,347)]
#    mlong=melt(mdf,id="pos")
#    max_prob=data.frame(pos=mdf$pos,total=rowSums(mdf[1:344]),founder=rep(founders[f],dim(mdf)[1]),chr=rep(c,dim(mdf)[1]))
#    coverage=rbind(coverage,max_prob)
#}
#cov_counts=coverage %>% group_by(pos,founder) %>% summarize(tot_cov=round(total,2))

#fwrite(cov_counts,sprintf('founder_coverage_counts_chr%s.txt',c),sep='\t',row.names=F,quote=F)

samples=dimnames(all_pr[[1]])[[1]]
options(warn=-1)
n_ind=length(samples)
xo_no=c()
all_breaks=c()
for (m in 1:n_ind){
  md = do.call(cbind,lapply(all_pr,function(x) x[m,]))
  #md <- all_pr_c[[1]][m,,] #pull out the probs for one individual
  #tm <- t(md) #transpose matrix
  mdf <- as.data.frame(md,row.names = rownames(md))
  mdf<-rownames_to_column(mdf,"marker")
  #names(mdf)=c("marker",founders)
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
    chr=c
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

fwrite(all_breaks,sprintf('stats/crossovers/Biogemma_breakpoints_chr%s.txt',c),row.names=F,quote=F,sep='\t')

xo=data.frame(sample=samples,chr=c,xo_no=xo_no)
fwrite(xo,sprintf('stats/crossovers/Biogemma_xo_number_chr%s.txt',c),row.names=F,quote=F,sep='\t')

#colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
#rownames(colorcodes)=colorcodes$founder

#founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#colorcodes=colorcodes[founders,]

#all_pr=c()
#plot_list=list()
#count=1
#for(c in 1:10){
#  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
#  pr=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
#  fsums=as.data.frame(lapply(pr,function(x) colSums(x)),stringsAsFactors=F)
#  fsums$marker=rownames(fsums)
#  fmelt=melt(fsums,'marker')
#  fmelt$chrom=c
#  fmelt$pos=pmap[match(fmelt$marker,pmap$marker),]$pos
#  if(c==1){
#    p=ggplot(fmelt,aes(x=pos/1e6,y=value,group=variable,color=variable)) + geom_line() + scale_color_manual(values=colorcodes$hex_color) + xlab("Position (Mb)") + ylab("Founder Coverage") + labs(subtitle=sprintf('Chromosome %.0f',c))
#    legend <- get_legend(
#      # create some space to the left of the legend
#      p + theme(legend.box.margin = margin(0, 0, 0, 12))
#    )
#    #p=p+guides(color=F,fill=F)
#  }else{
#    p=ggplot(fmelt,aes(x=pos/1e6,y=value,group=variable,color=variable)) + geom_line() + scale_color_manual(values=colorcodes$hex_color) + xlab("Position (Mb)") + ylab("Founder Coverage") + labs(subtitle=sprintf('Chromosome %.0f',c))
#  }
#  plot_list[[count]]=p
#  count=count+1
#}

#allp=plot_grid(plotlist=plot_list,nrow=10)
#png('stats/coverage/coverage.png',width=2000,height=1000)
#print(allp)
#dev.off()

#pdf('stats/coverage/coverage.pdf',width=14)
#for(i in 1:length(plot_list)){
#  print(plot_list[[i]])
#}
#dev.off()
