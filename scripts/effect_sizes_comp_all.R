#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
ft_days=c('male_flowering_days','female_flowering_days')
qtl_overlap=qtl_overlap[!(qtl_overlap$Phenotype %in% ft_days),]
rownames(qtl_overlap)=seq(1,nrow(qtl_overlap))
#qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)

ses=readRDS('GridLMM/effect_sizes/All_QTL_ES.rds')

make_plot<-function(data,name,focal,grouptype){
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[data$founder,]
  leg<-ggplot(data,aes(x=variable_f,y=f_value,fill=variable_f)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=colorcodes[levels(data$variable_f),]$hex_color,labels=levels(data$variable_f)) +
  guides(fill=guide_legend(title="Founder")) +
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6))


  legend <- get_legend(
    leg + theme(legend.box.margin = margin(0, 0, 0, 20))
  )

  a<-ggplot(data,aes(x=variable_f,y=f_value,color=variable_f)) +
  geom_point() + scale_color_manual(values=colorcodes[levels(data$variable_f),]$hex_color,labels=levels(data$variable_f)) +
  geom_errorbar(aes(ymin=f_value-(2*f_se),ymax=f_value+(2*f_se)),width=.2,position=position_dodge()) +
  geom_text(aes(label=allele),vjust=-0.3,hjust=-0.3,color="black",size=3.5) +
  ylab("Effect Size") + xlab("Founder") +
  ggtitle(sprintf('%s Effect Sizes of %s, %s Focal SNP',grouptype,name,focal)) +
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5),axis.title.y=element_text(size=6),axis.title.x=element_text(size=8),plot.title=element_text(size=12)) +
  guides(color=F)
  sdata = data %>% group_by(allele) %>% summarize(s_value=mean(s_value),s_se=mean(s_se))

  b<-ggplot(sdata,aes(x=factor(allele,levels=c(0,1)),y=s_value)) +
  geom_point() + geom_errorbar(aes(ymin=s_value-(2*s_se),ymax=s_value+(2*s_se))) +
  ylab("SNP Effect Size") +
  xlab("Allele") +
  #ggtitle(sprintf('SNP Effect Sizes of %s',name)) +
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5),
  axis.title.x=element_text(size=8),axis.title.y=element_text(size=6))

  hdata = data %>% group_by(hapgrp) %>% summarize(h_value=mean(h_value),h_se=mean(h_se))
  hdata$variable_h=factor(paste0("HAPGRP_",hdata$hapgrp),levels=paste0("HAPGRP_",hdata$hapgrp))
  c<-ggplot(hdata,aes(x=variable_h,y=h_value)) + geom_point() +
  geom_errorbar(aes(ymin=h_value-(2*h_se),ymax=h_value+(2*h_se))) +
  ylab("Haplotype Effect Size") +
  xlab("Haplotype Groups") +
  theme(axis.text.x=element_text(size=5),
  axis.text.y=element_text(size=5),
  axis.title.y=element_text(size=4),
  axis.title.x=element_text(size=8))

  c2=list()
  count=1
  h=max(data$hapgrp)
  for(i in seq(1,h)){
    sub=data[data$hapgrp==i,]
    c2[[count]]<- ggplot(sub,aes(x="",y=h_perc,fill=variable_f))+
    geom_bar(stat="identity",width=1) +
    coord_polar("y",start=0) +
    scale_fill_manual(values=colorcodes[sub$variable_f,]$hex_color,labels=levels(sub$variable_f)) +
    guides(fill=F) +
    labs(caption=paste0('HAPGRP_',i)) +
    theme(axis.title.x=element_blank(),
    axis.text.y=element_blank(),axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),plot.caption=element_text(size=4))

    count=count+1
  }

  pies=plot_grid(plotlist=c2,nrow=1)
  p1=plot_grid(pies,c,ncol=1,rel_heights=c(3,12),rel_widths=c(1,1.4))
  p2=plot_grid(b,p1,nrow=1,rel_widths=c(3,h))
  p3=plot_grid(a,p2,nrow=2,rel_heights=c(8,8))
  p4=plot_grid(p3,legend,ncol=2,rel_widths=c(4,1.2))
  return(p4)
}
s_f_h_data=list()
s_f_h_plots=list()
acount=1
l=length(ses)
for(i in seq(1,l,3)){
  name=ses[[i]]$qtl_id
  chr=unique(qtl[qtl$pheno_env_id==name,]$Chromosome)
  ssnp=ses[[i]]$keptsnp
  smelt=ses[[i]]$SE
  fmelt=ses[[i+1]]$SE
  if(nrow(fmelt)<16){
    rown=rownames(fmelt)
    naf=founders[!founders %in% rownames(fmelt)]
    for(f in naf){
      line=rep(NA,5)
      #rownames(line)=f
      fmelt=rbind(fmelt,line)
    }
    rownames(fmelt)=c(rown,naf)
    fmelt=fmelt[founders,]
  }
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')
  hmelt=ses[[i+2]]$SE
  hsnp=ses[[i+2]]$keptsnp
  h=ses[[i+2]]$hapgrp
  allhaps=paste0('HAPGRP_',seq(1,h))
  if(nrow(hmelt)<h){
    rown=rownames(hmelt)
    naf=allhaps[!allhaps %in% rownames(hmelt)]
    for(f in naf){
      line=rep(NA,5)
      #rownames(line)=f
      hmelt=rbind(hmelt,line)
    }
    rownames(hmelt)=c(rown,naf)
    hmelt=hmelt[allhaps,]
  }

  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=ses[[i]]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps
  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)

  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  names(smelt)=c('value','se','tvalue','allele')

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue

  grouptype=qtl_overlap[qtl_overlap$pheno_env_id==name,]$label
  focal=ses[[i]]$focal
  s_f_h_data[[acount]]=list(data=fmelt,qtl_id=name,focal=focal,grouptype=grouptype)
  s_f_h_plots[[acount]]=make_plot(fmelt,name,focal,grouptype)
  acount=acount+1
}

saveRDS(s_f_h_data,'GridLMM/effect_sizes/Focal_QTL_Effect_Sizes.rds')


pdf('GridLMM/effect_sizes/Focal_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(s_f_h_plots)){
  print(s_f_h_plots[[i]])
}
dev.off()
