#!/usr/bin/env Rscript

library('data.table')
library('tidyverse')

qtl=fread('Biogemma_QTL_update.csv',data.table=F)
base_list=c(8,7,7,8,6,7,7,8,7,7)

alt_right_bound_bp=c()
alt_right_bound_snp=c()

for (i in seq(1,dim(qtl)[1])) {
  if (qtl[i,]$Method == "Founder_probs") {
    print('Founder_prob')
    chr=qtl[i,]$Chromosome
    m=qtl[i,]$right_bound_snp
    fdropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%.0f_dropped_markers_genoprobs.rds',chr))
    dropped_markers=fdropped[[which(unlist(lapply(fdropped,function(x) x$marker==m)))]]$linked
    pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
    sub=pmap[pmap$marker %in% dropped_markers,]
    rownames(sub)=seq(1,dim(sub)[1])
    right=which.max(sub$pos)
    right_bound=sub[right,]$pos
    right_snp=sub[right,]$marker
    #print(c(i,length(right_bound)))
    alt_right_bound_bp=c(alt_right_bound_bp,right_bound)
    alt_right_bound_snp=c(alt_right_bound_snp,right_snp)
  }
  else if (qtl[i,]$Method == "Haplotype_probs") {
    print('Haplotype_prob')
    chr=qtl[i,]$Chromosome
    m=qtl[i,]$right_bound_snp
    for (h in seq(base_list[chr],16)) {
      hdropped=readRDS(sprintf('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%.0f_haplogroup%.0f_dropped_markers.rds',chr,h))
      if (length(which(unlist(lapply(hdropped,function(x) x$marker==m)))) != 0) {
        print(h)
        dropped_markers=hdropped[[which(unlist(lapply(hdropped,function(x) x$marker==m)))]]$linked
        if (!is.null(dropped_markers)) {
          print('not null')
          pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
          sub=pmap[pmap$marker %in% dropped_markers,]
          rownames(sub)=seq(1,dim(sub)[1])
          right=which.max(sub$pos)
          right_bound=sub[right,]$pos
          right_snp=sub[right,]$marker
          #print(c(i,length(right_bound)))
          alt_right_bound_bp=c(alt_right_bound_bp,right_bound)
          alt_right_bound_snp=c(alt_right_bound_snp,right_snp)
        }
        else {
          print('Null')
          alt_right_bound_bp=c(alt_right_bound_bp,qtl[i,]$right_bound_bp)
          alt_right_bound_snp=c(alt_right_bound_snp,qtl[i,]$right_bound_snp)
        }
      }
    }
  }
  else {
    print('600K_SNP')
    #print(c(i,length(qtl[i,]$right_bound_bp)))
    alt_right_bound_bp=c(alt_right_bound_bp,qtl[i,]$right_bound_bp)
    alt_right_bound_snp=c(alt_right_bound_snp,qtl[i,]$right_bound_snp)
  }
  print(c(i,length(alt_right_bound_bp)))
}

qtl$alt_right_bound_bp=alt_right_bound_bp
qtl$alt_right_bound_snp=alt_right_bound_snp
qtl$alt_size=qtl$alt_right_bound_bp - qtl$left_bound_bp

chrom_lengths=c(306972432,244416191,235604081,246907405,223714093,173392688,182370476,181098567,159740209,150903715)

cum_left=c()
cum_right=c()
alt_cum_right=c()

for(i in 1:dim(qtl)[1]){
  print(i)
  left_pos=qtl[i,]$left_bound_bp
  right_pos=qtl[i,]$right_bound_bp
  alt_right_pos=qtl[i,]$alt_right_bound_bp
  total=0
  chr=qtl[i,]$Chromosome
  if(chr==1){
    cum_left=c(cum_left,left_pos)
    cum_right=c(cum_right,right_pos)
    alt_cum_right=c(alt_cum_right,alt_right_pos)
  }
  else{
    for(c in 1:(qtl[i,]$Chromosome-1)){
      total=total+chrom_lengths[c]
    }
    cum_left_pos=total+left_pos
    cum_left=c(cum_left,cum_left_pos)
    cum_right_pos=total+right_pos
    cum_right=c(cum_right,cum_right_pos)
    alt_cum_right_pos=total+alt_right_pos
    alt_cum_right=c(alt_cum_right,alt_cum_right_pos)
  }
}
qtl$cum_left_bound_bp=cum_left
qtl$cum_right_bound_bp=cum_right
qtl$alt_cum_right_bound_bp=alt_cum_right

fwrite(qtl,'Biogemma_QTL_update.csv',row.names=F,quote=F,sep=',')
