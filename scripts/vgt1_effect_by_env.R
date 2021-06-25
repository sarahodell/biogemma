library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

envs=c("ALL","BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD","SZEGED_2017_OPT","STPAUL_2017_WD")
names=paste0("male_flowering_d6_",envs,"_qDTA8")
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
has_mite=list("B73_inra"=F,"A632_usa"=T,"CO255_inra"=T,"FV252_inra"=T,"OH43_inra"=F,
           "A654_inra"=T,"FV2_inra"=T,"C103_inra"=T,"EP1_inra"=T,"D105_inra"=T,
           "W117_inra"=T,"B96"=F,"DK63"=T,"F492"=T,"ND245"=T,"VA85"=F)

plot_list=list()
count=1
chr=8
ses=readRDS('GridLMM/effect_sizes/All_QTL_ES.rds')

for(name in names){
  n=which(unlist(lapply(ses,function(x) x$qtl_id == name & x$focal=="Founder_probs")))
  #i=n[1]
  ssnp=ses[[n[1]]]$keptsnp
  smelt=ses[[n[1]]]$SE
  fmelt=ses[[n[2]]]$SE
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
  hmelt=ses[[n[3]]]$SE
  hsnp=ses[[n[3]]]$keptsnp
  h=ses[[n[3]]]$hapgrp
  allhaps=paste0('HAPGRP_',seq(1,h))
  if(nrow(hmelt)<h){
    rown=rownames(hmelt)
    naf=allhaps[!allhaps %in% rownames(hmelt)]
    for(f in naf){
      line=rep(NA,5)
      hmelt=rbind(hmelt,line)
    }
    rownames(hmelt)=c(rown,naf)
    hmelt=hmelt[allhaps,]
  }
  qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=ses[[n[1]]]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%.0f_genoprobs.rds',chr))
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
  focal=ses[[n[1]]]$focal
  fmelt$mite=unlist(unname(sapply(seq(1,16),function(x) has_mite[rownames(fmelt)[x]])))
  founderlabels=c("B73","OH43","A632","CO255","F252","A654","DK63","F2","EP1","D105","W117","F492","C103","B96","ND245","VA85")
  hap_means=fmelt %>% group_by(hapgrp) %>% summarize(mean=mean(f_value,na.rm=T),n=length(f_value))
  fmelt$hap_value=hap_means[fmelt$hapgrp,]$mean
  fmelt=fmelt[order(fmelt$mite,fmelt$hap_value),]
  fmelt$hap_allele=paste0(fmelt$allele,'_',fmelt$hapgrp)
  fmelt$variable_f=factor(fmelt$founder,levels=fmelt$founder)
  dup=which(duplicated(fmelt$hap_allele))
  i=seq(1,16,2)
  f=rownames(fmelt)[i]
  uniq=length(unique(fmelt$hap_allele))
  for(d in dup){
    #if(rownames(fmelt)[d] %in% f){
    i[i>=d]=i[i>=d]+1
    i=i[i<=16]
    #f=rownames(fmelt)[i]
      #m=which(fmelt$hap_allele==fmelt$hap_allele[d])
  }



  df <- data.frame(group = fmelt$hap_allele[i],
                   variable_f = factor(fmelt$variable_f[i],levels=levels(fmelt$variable_f)),
                   Y = seq(1,length(i)))
  df$xmin=as.numeric(fmelt$variable_f[i])-0.5
  df$xmax=as.numeric(fmelt$variable_f[i])+0.5
  for(d in dup){
    if((d-1) %in% i){
      df[df$variable_f %in% fmelt$variable_f[d-1],]$xmax=df[df$variable_f %in% fmelt$variable_f[d-1],]$xmax+1
    }
  }


  hjust=rep(0.5,16)
  hjust[(dup-1)]=-2
  hjust=hjust[-dup]
  df2=data.frame(group=unique(fmelt$hap_allele),
    variable_f=fmelt$variable_f[-dup],
    hjust=hjust,
    vjust=25)
  scale=seq(-60,60,10)
  ymin=min(fmelt$f_value-(2*fmelt$f_se),na.rm=T)
  ymax=max(fmelt$f_value+(2*fmelt$f_se),na.rm=T)
  lower=scale[(scale-ymin)>(-10) & (scale-ymin<0)]
  upper=scale[(scale-ymax)<(10) & (scale-ymax>0)]
  scale=seq(lower,upper,10)
  if(name=="male_flowering_d6_ALL_qDTA8"){
    a<-ggplot() +
     geom_point(data=fmelt,aes(x=variable_f,y=f_value,color=allele,group=hap_allele),position = position_dodge(width=1),size=10) +
     geom_errorbar(data=fmelt,aes(x=variable_f,ymin=f_value-(2*f_se),ymax=f_value+(2*f_se),color=allele),position = position_dodge(width=1),size=4,width=.5) +
     geom_rect(data=df,aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),fill="grey60",alpha=0.4) +
     geom_text(data=df2,aes(x=variable_f,y=38,hjust=hjust),label=LETTERS[1:uniq],color="black",size=14)+
     geom_hline(yintercept=0,color="black") +
     ylab("DTA Effect Size (ggd)") + xlab("Founder") +
     labs(color="Allele") +
     scale_x_discrete(labels=c("B73","OH43","VA85","B96","EP1","ND245","D105","A654","DK63","F2","CO255","W117","F492","C103_inra","A632","F252")) +
     scale_y_continuous(limits=c(lower,upper),breaks=scale)+
     theme_classic() +
     #ggtitle(sprintf("Founder Effect Sizes %s",name)) +
     theme(axis.text.y=element_text(size=24),
     axis.title.y=element_text(size=30),
     axis.title.x=element_text(size=30),
     legend.title = element_text(size = 40),
     legend.position = c(0.9,0.1),
     legend.text = element_text(size = 36),
     panel.grid.major.y=element_line(color="black"),
      axis.text.x=element_text(size=24,vjust=0.5)) +
     theme(axis.text.x=element_text(face=ifelse(levels(fmelt$variable_f) %in% c("EP1_inra","D105_inra","A654_inra","DK63","FV2_inra","CO255_inra","W117_inra","F492","A632_usa","FV252_inra","ND245","C103_inra"),"bold.italic","plain")))
  }else{
    a<-ggplot() +
     geom_point(data=fmelt,aes(x=variable_f,y=f_value,color=allele,group=hap_allele),position = position_dodge(width=1),size=6) +
     geom_errorbar(data=fmelt,aes(x=variable_f,ymin=f_value-(2*f_se),ymax=f_value+(2*f_se),color=allele),position = position_dodge(width=1),size=2,width=.5) +
     geom_rect(data=df,aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),fill="grey60",alpha=0.4) +
     #geom_text(data=df2,aes(x=variable_f,y=35,hjust=hjust),label=LETTERS[1:uniq],color="black",size=6)+
     geom_hline(yintercept=0,color="black") +
     ylab("DTA Effect Size (ggd)") + xlab("Founder") +
      labs(color="Allele") +
      scale_y_continuous(limits=c(lower,upper),breaks=scale)+
     theme_classic() +
     ggtitle(sprintf("Founder Effect Sizes %s",name)) +
     theme(axis.text.y=element_text(size=12),
     axis.title.y=element_text(size=12),
     axis.title.x=element_text(size=12),
     legend.title = element_text(size = 14),
     legend.text = element_text(size = 12),
     panel.grid.major.y=element_line(color="black"),
      axis.text.x=element_text(size=10,angle=-45,vjust=0.5)) +
     theme(axis.text.x=element_text(face=ifelse(levels(fmelt$variable_f) %in% c("EP1_inra","D105_inra","A654_inra","DK63","FV2_inra","CO255_inra","W117_inra","F492","A632_usa","FV252_inra","ND245","C103_inra"),"bold.italic","plain")))
  }
   plot_list[[count]]=a
   count=count+1
}

png('GridLMM/effect_sizes/vgt1_Figure4.png',width=2000,height=1250)
print(plot_list[[1]])
dev.off()

pdf('GridLMM/effect_sizes/Methods_Supplemental3.pdf',width=10,height=8)
for(l in 2:length(plot_list)){
  print(plot_list[[l]])
}
dev.off()
