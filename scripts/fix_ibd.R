#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')
founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
hap_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
ibd=fread(sprintf('ibd_segments/bg%s_wgs_ibdsegments_010319.txt',c),data.table=F)

pmap=fread(sprintf('qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)
fgeno=fread(sprintf('../Biogemma_foundergenos/Founder_genos_chr%s_121718.csv',c),data.table=F)
pmap$pos=round(pmap$pos)
rownames(fgeno)=founders
fgeno=fgeno[hap_founders,c('ind',pmap$marker)]
markers=unlist(colnames(fgeno))
new_ibd=c()

for(i in seq(1,dim(ibd)[1])){
      line=ibd[i,]
      leftind=line$left_index+1
      rightind=leftind+line$n_mar-2
      leftm=pmap[leftind,]$marker
      rightm=pmap[rightind,]$marker
      leftp=pmap[leftind,]$pos
      rightp=pmap[rightind,]$pos
      intlen=rightp-leftp

      # Find number of mismatches
      left=which(markers==leftm)
      right=which(markers==rightm)
      bin=fgeno[c(line$strain1,line$strain2),left:right]
      mm=sum(apply(bin,MARGIN=2,FUN=function(x) length(unique(x)))>1)
      new_line=c(line$strain1,line$strain2,line$chr,leftm,rightm,leftp,rightp,leftind,rightind,intlen,line$n_mar,mm,line$lod)
      new_ibd=rbind(new_ibd,new_line)

}

new_ibd=as.data.frame(new_ibd,stringsAsFactors=F)
names(new_ibd)=c('strain1','strain2','chr','left_marker','right_marker','left_pos','right_pos','left_index','right_index','int_length','n_mar','n_mismatch','lod')

print(sum(new_ibd$n_mismatch>ibd$n_mismatch))

fwrite(new_ibd,sprintf('ibd_segments/bg%s_ibd_segments_adjusted.txt',c),row.names=F,quote=F,sep='\t')