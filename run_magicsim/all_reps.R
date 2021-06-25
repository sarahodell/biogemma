#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')


all_sig=c()
#sig=c()
for(r in 1:100){
  total_tests=0
  tmpsig=c()
  for(c in 1:10){
    pchi=fread(sprintf('selection/bg%.0f_founder_chisq_rep%.0f_results.txt',c,r),data.table=F)
    pchi$chr=c
    tmpsig=rbind(tmpsig,pchi)
    total_tests=total_tests+nrow(pchi)
  }
  total_tests=total_tests*16
  bonf=-log10(0.05 / total_tests)
  #s=tmpsig$p_chi[-log10(tmpsig$p_chi)>=bonf]
  t=tmpsig[-log10(tmpsig$p_chi)>=bonf,]
  if(nrow(t)!=0){
    #t$chr=c
    t$rep=r
    all_sig=rbind(all_sig,t)
  }
  #sig=rbind(sig,c(r,length(s)))
}
sig=as.data.frame(sig,stringAsFactors=F)
names(sig)=c('rep','sig_count')

fwrite(sig,'selection/total_sim_chi_sq_peaks.txt',row.names=F,quote=F,sep='\t')

png('selection/chi_peak_sim_dist.png')
print(ggplot(sig,aes(x=sig_count)) + geom_histogram(bins=20) + xlab("Number of Chi Squared Peaks") + ylab("Frequency"))
dev.off()

# Get peak size, percentage of genome

segments=c()
#baselist=c(8,7,7,8,6,7,7,8,7,7)
for(i in 1:nrow(all_sig)){
  line=all_sig[i,]
    #h=line$hapgrp
  left_bound=line$pos
  m=line$marker
  chr=line$chr
  rep=line$rep
  fdropped=readRDS(sprintf('qtl2_files/dropped/bg%.0f_rep%.0f_dropped_markers.rds',chr,rep))
  index=which(unlist(lapply(fdropped,function(x) x$marker==m)))
  dropped_markers=fdropped[[index]]$linked
  if(!is.null(dropped_markers)){
    pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
    sub=pmap[pmap$marker %in% dropped_markers,]
    rownames(sub)=seq(1,dim(sub)[1])
    right=which.max(sub$pos)
    right_snp=sub[right,]$marker
    right_bound=sub[right,]$pos
  }else{
    right_bound=pmap[(which(pmap$marker==m)+1),]$pos - 1
    right_snp=all_sig[(i+1),]$marker
  }
  seg=c(rep,chr,left_bound,right_bound,paste0(m,'..',right_snp))
  segments=rbind(segments,seg)
}


segments=as.data.frame(segments,stringsAsFactors=F)
rownames(segments)=seq(1,nrow(segments))
names(segments)=c('rep','chr','start','end','markers')
segments$rep=as.numeric(segments$rep)
segments$chr=as.numeric(segments$chr)
segments$start=as.numeric(segments$start)
segments$end=as.numeric(segments$end)
fwrite(segments,'selection/chi_peak_seqments.bed',row.names=F,quote=F,sep='\t',col.names=F)


segments$size=segments$end-segments$start
total=0
for(c in 1:10){
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  total=total+max(pmap$pos)
}
perc= segments %>% group_by(rep) %>% summarize(chi_perc=sum(size)/total)
