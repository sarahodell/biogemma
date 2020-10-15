#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')

total_tests=125764
bonf=-log10(0.05 / total_tests)
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
chipeaks=fread(sprintf('selection/haplotype_probs/bg%s_haplotype_chisq_results.txt',chr),data.table=F)
chipeaks$sig=-log10(chipeaks$p_chi)>=bonf

segments=c()
baselist=c(8,7,7,8,6,7,7,8,7,7)
for(i in 1:nrow(chipeaks)){
  line=chipeaks[i,]
  if(line$sig){
    h=line$hapgrp
    left_bound=line$pos
    m=line$marker
    hdropped=readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/dropped/bg%s_haplogroup%.0f_dropped_markers.rds',chr,h))
    dropped_markers=hdropped[[which(unlist(lapply(hdropped,function(x) x$marker==m)))]]$linked
    if (!is.null(dropped_markers)) {
      sub=pmap[pmap$marker %in% dropped_markers,]
      rownames(sub)=seq(1,dim(sub)[1])
      right=which.max(sub$pos)
      right_snp=sub[right,]$marker
      right_bound=sub[right,]$pos
    }
    else {
      right_bound=pmap[(which(pmap$marker==m)+1),]$pos - 1
      right_snp=chipeaks[(i+1),]$marker
    }
    seg=c(chr,left_bound,right_bound,paste0(m,'..',right_snp),paste0('hapgrp',h))
    segments=rbind(segments,seg)
  }
}

segments=as.data.frame(segments)
rownames=seq(1,dim(segments)[1])
names(segments)=c('chr','start','end','markers','hapgrp')
fwrite(segments,sprintf('selection/haplotype_probs/bg%s_chi_peak_seqments.bed',chr),row.names=F,quote=F,sep='\t',col.names=F)
