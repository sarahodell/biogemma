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
qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)

#f_ses=readRDS('GridLMM/effect_sizes/Founder_prob_QTL_SEs.rds')
h_ses=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
#s_ses=readRDS('GridLMM/effect_sizes/600K_SNP_QTL_SEs.rds')

#s_notf_ses=readRDS('GridLMM/effect_sizes/SNP_notFounder_prob_QTL_SEs.rds')
h_notf_ses=readRDS('GridLMM/effect_sizes/Haplotype_notFounder_prob_QTL_SEs.rds')
#f_nots_ses=readRDS('GridLMM/effect_sizes/Founder_notSNP_prob_QTL_SEs.rds')

f_nots=fread('GridLMM/effect_sizes/Founder_notSNP_highest_peaks.txt',data.table=F)
s_notf=fread('GridLMM/effect_sizes/SNP_notFounder_highest_peaks.txt',data.table=F)

f_noth=fread('GridLMM/effect_sizes/Founder_notHaplotype_highest_peaks.txt',data.table=F)
h_notf=fread('GridLMM/effect_sizes/Haplotype_notFounder_highest_peaks.txt',data.table=F)
# Founder vs. Haplotype for shared QTL
ft_days=c("female_flowering_days","male_flowering_days")
qtl_overlap=qtl_overlap[!(qtl_overlap$Phenotype %in% ft_days),]
q1=qtl_overlap[qtl_overlap$F_and_H !=0,]
m1=c('Founder_probs','Haplotype_probs','600K_SNP')

f_h = qtl[qtl$pheno_env %in% q1$pxe,]
f_h$pheno_env_id = paste0(f_h$pheno_env, '_',f_h$ID)
counts=f_h %>% group_by(pheno_env_id) %>% count
keep=counts[counts$n==2,]$pheno_env_id

f_h = f_h[f_h$pheno_env_id %in% keep,]
rownames(f_h)=seq(1,nrow(f_h))

f_h_ids=unique(f_h$pheno_env_id)

f_h_plots=list()
count=0

for(q in f_h_ids){
  sub=f_h[f_h$pheno_env_id == q,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method==m1[1],]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fmelt=as.data.frame(f_ses[[which(sapply(f_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
  fmelt$variable=founders
  rownames(fmelt)=founders
  #fmelt$method="Founder_probs"
  #fmelt$hapgrp=factor(seq(1,16),levels=seq(1,16))
  names(fmelt)=c('f_value','f_se','f_tvalue','variable')
  #fmelt$variable_f=factor(fmelt$variable,levels=fmelt[order(fmelt$hapgrp),]$variable)

  #Haplotype effect sizes
  line=sub[sub$Method==m1[2],]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=hap_table[hap_table$SNP == line$highest_SNP, ]$HAPGRP
  hmelt=as.data.frame(h_ses[[which(sapply(h_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==line$highest_SNP,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$hapgrp=factor(ibd_seg,levels=seq(1,h))
  names(hmelt)=c('h_value','h_se','h_tvalue','variable','hapgrp')
  hmelt$variable_f=factor(hmelt$variable,levels=hmelt[order(hmelt$hapgrp),]$variable)
  #a<-ggplot(hmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))
  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_f = hmelt$variable_f

  line=f_nots[f_nots$pheno_env_id==q,]
  smelt=as.data.frame(f_nots_ses[[which(sapply(f_nots_ses, function(x) x$qtl_id==q))]]$SE,stringsAsFactors=F)
  ssnp=line$highest_SNP
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f_121718.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[founders,ssnp]=="A",2,1)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  #smelt$method="600K_SNP"
  alleles=ifelse(geno[founders,ssnp]=="A",0,1)
  smelt$allele=alleles
  #smelt$variable=as.character(smelt$variable)
  names(smelt)=c('value','se','tvalue','variable','allele')
  #allmelt=rbind(hmelt,fmelt,smelt)

  fmelt$allele=smelt$allele
  fmelt$s_value=sapply(seq(1,16), function(x) ifelse(!is.na(smelt$value[x]),smelt$value[x],0))
  #count=count+1
  #f_h_plots[[count]]=ggplot(allmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))
  count=count+1
  f_h_plots[[count]]=list(values=fmelt,id=q)
}

saveRDS(f_h_plots,'GridLMM/effect_sizes/F_and_H_all_effect_sizes.rds')

f_h_data=readRDS('GridLMM/effect_sizes/F_and_H_all_effect_sizes.rds')
l=length(f_h_data)
f_h_plots=list()
count=0
for(n in seq(1,l)){
  #print(n)
  data=as.data.frame(f_h_data[[n]]$values,stringsAsFactors=F)
  name=f_h_data[[n]]$id
  #print(name)
  #ylim=max(data$f_value + f_se)
  data$variable_f=factor(data$variable,levels=data[order(data$hapgrp),]$variable_f)
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[data$variable,]
  a<-ggplot(data,aes(x=variable,y=f_value)) +
   geom_bar(aes(fill=variable),stat="identity") +
    scale_fill_manual(values=colorcodes[data$variable,]$hex_color,labels=data$variable) +
     geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge()) +
      facet_grid(.~allele) + ylab("Effect Size") + xlab("Founder") +
       ggtitle(sprintf('Founder BLUP Effect Sizes of %s',name)) +
        theme(axis.text.x=element_text(size=8))
  legend <- get_legend(
  # create some space to the left of the legend
    a + theme(legend.box.margin = margin(0, 0, 0, 20))
  )


  data0=data[data$allele==0,]
  a0<-ggplot(data0,aes(x=variable,y=f_value)) +
   geom_bar(aes(fill=variable),stat="identity") +
    scale_fill_manual(values=colorcodes[data0$variable,]$hex_color,labels=data0$variable) +
     geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge()) +
       ylab("Founder Effect Size (ggd)") + xlab("Founder")  + theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) +
       labs(subtitle="Allele 0 (Reference) of 600K SNP") + guides(fill=F) #+  theme_classic()
  #a0

  data1=data[data$allele==1,]
  a1<-ggplot(data1,aes(x=variable,y=f_value)) +
   geom_bar(aes(fill=variable),stat="identity") +
    scale_fill_manual(values=colorcodes[data1$variable,]$hex_color,labels=data1$variable) +
     geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge()) + ylab("Founder Effect Size") +
      xlab("Founder") + theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank()) +
       labs(subtitle="Allele 1 (Alternate) of 600K SNP") + guides(fill=F)

  hdata = data %>% group_by(hapgrp) %>% summarize(h_value=mean(h_value),h_se=mean(h_se))
  hdata=as.data.frame(hdata,stringsAsFactors=F )
  hdata$hapgrp=as.numeric(hdata$hapgrp)
  h=max(hdata$hapgrp)
  founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
  founder_name=c(rep(founders,h))
  haplos=c()
  for(i in seq(1,h)){haplos=c(haplos,rep(i,16))}
  founder_prop=c()
  for(j in seq(1,h)){
    sub=data[data$hapgrp==j,]
    f=unique(sub$variable)
    no=length(f)
    for(k in founders){
      if(k %in% f){
        founder_prop=c(founder_prop,1/no)
      }
      else{
        founder_prop=c(founder_prop,0)
      }
    }
  }
  fprop=data.frame(founder=founder_name,hapgrp=haplos,founder_prop=founder_prop)

  h_value=c()
  h_se=c()
  for(o in seq(1,h)){
    v=hdata[hdata$hapgrp==o,]$h_value
    se=hdata[hdata$hapgrp==o,]$h_se

    h_value=c(h_value,rep(v,16))
    h_se=c(h_se,rep(se,16))
  }
  fprop$h_value=h_value
  fprop$h_se=h_se
  fprop$h_value_adj=fprop$founder_prop * fprop$h_value
  b<-ggplot(fprop,aes(x=factor(hapgrp,levels=seq(1,h)),y=h_value_adj,fill=founder)) +
   geom_bar(position="stack",stat="identity")+ scale_fill_manual(values=colorcodes[data$variable,]$hex_color,labels=data$variable)+
    ylab("Haplotype Effect Size") +
    geom_errorbar(aes(ymin=h_value-h_se,ymax=h_value+h_se),width=.02,position=position_dodge()) +
    xlab("Haplotype Groups") + ggtitle(sprintf('Haplotype BLUP Effect Sizes of %s',name)) +
     theme(axis.text.x=element_text(size=8)) + guides(fill=F)

  p1=plot_grid(a0,a1,ncol=2,rel_widths=c(9,6)) + xlab("Founders") + labs(title=sprintf("Founder BLUP Effect Sizes of %s",name))
  p2=plot_grid(p1,b,nrow=2)
  p3=plot_grid(p2,legend,ncol=2,rel_widths = c(4,.4))
  p4=p3 + labs(title=sprintf("BLUP Effect Sizes of %s",name))
  count=count+1
  f_h_plots[[count]]=p4
}

pdf('GridLMM/effect_sizes/F_and_H_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(f_h_plots)){
  print(f_h_plots[[i]])
}
dev.off()
