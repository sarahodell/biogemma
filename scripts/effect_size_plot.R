library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
ft_days=c("female_flowering_days","male_flowering_days")
qtl=qtl[!(qtl$Phenotype %in% ft_days),]
qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
qtl_overlap=qtl_overlap[!(qtl_overlap$Phenotype %in% ft_days),]

#qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)

f_ses=readRDS('GridLMM/effect_sizes/founder_ES/Founder_prob_QTL_SEs.rds')
h_ses=readRDS('GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')
s_ses=readRDS('GridLMM/effect_sizes/600K_SNP_QTL_SEs.rds')
# Founder vs. Haplotype for shared QTL
q1=qtl_overlap[qtl_overlap$label=="S_F_H",]
m1=c('Founder_probs','Haplotype_probs',"600K_SNP")

#s_f_h = qtl[qtl$pheno_env %in% q1$pxe,]
#s_f_h$pheno_env_id = paste0(s_f_h$pheno_env, '_',s_f_h$ID)
#counts=s_f_h %>% group_by(pheno_env_id) %>% count
#keep=counts[counts$n==3,]$pheno_env_id
#s_f_h = s_f_h[s_f_h$pheno_env_id %in% keep,]
#rownames(s_f_h)=seq(1,nrow(s_f_h))

s_f_h_ids=unique(q1$pheno_env_id)

s_f_h_plots=list()
count=0

for(q in s_f_h_ids){
  sub=qtl[qtl$pheno_env_id == q,]
  rownames(sub)=seq(1,nrow(sub))
  #Founder effect sizes
  line=sub[sub$Method==m1[1],]
  pheno=line$Phenotype
  env=line$Environment
  chr=line$Chromosome
  fdata=f_ses[[which(sapply(f_ses, function(x) x$qtl_id==q))]]
  fmelt=as.data.frame(fdata$SE,stringsAsFactors=F)
  fsnp=fdata$snp
  #fmelt$variable=founders
  rownames(fmelt)=founders
  #fmelt$method="Founder_probs"
  #fmelt$hapgrp=factor(seq(1,16),levels=seq(1,16))
  names(fmelt)=c('f_value','f_se','f_tvalue','founder','variable_f')
  #fmelt$variable_f=factor(fmelt$variable,levels=fmelt[order(fmelt$hapgrp),]$variable)

  #Haplotype effect sizes
  line=sub[sub$Method==m1[2],]
  hap_table=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=hap_table[hap_table$SNP == line$highest_SNP, ]$HAPGRP
  hdata=h_ses[[which(sapply(h_ses, function(x) x$qtl_id==q))]]
  hsnp=hdata$snp
  hmelt=as.data.frame(hdata$SE,stringsAsFactors=F)
  ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==line$highest_SNP,]$pos
  ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))
  hreps=round(colSums(fprobs[[1]][,,hsnp]))
  names(hreps)=founders
  hmelt$hapgrp=seq(1,h)
  hmelt=hmelt[ibd_seg,]
  rownames(hmelt)=seq(1,nrow(hmelt))
  hmelt$variable=founders
  hmelt$reps=hreps

  
  names(hmelt)=c('h_value','h_se','h_tvalue','hapgrp','variable_f','variable','reps')
  hmelt = hmelt %>% group_by(hapgrp) %>% mutate(hap_total=sum(reps))
  hmelt$perc=round(hmelt$reps/hmelt$hap_total * 100)
  #hmelt$variable_f=factor(hmelt$variable,levels=hmelt[order(hmelt$hapgrp),]$variable)
  #a<-ggplot(hmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))
  fmelt$h_value=hmelt$h_value
  fmelt$h_se=hmelt$h_se
  fmelt$h_tvalue=hmelt$h_tvalue
  fmelt$hapgrp=hmelt$hapgrp
  fmelt$variable_h = hmelt$variable_f
  fmelt$h_perc=hmelt$perc

  line=sub[sub$Method==m1[3],]
  sdata=s_ses[[which(sapply(s_ses, function(x) x$qtl_id==q))]]
  smelt=as.data.frame(sdata$SE,stringsAsFactors=F)
  ssnp=sdata$snp
  geno=fread(sprintf('genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',chr),data.table=F)
  rownames(geno)=geno$ind
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",1,2)
  smelt=smelt[alleles,]
  rownames(smelt)=seq(1,nrow(smelt))
  smelt$variable=founders
  #smelt$method="600K_SNP"
  alleles=ifelse(geno[match(founders,geno$ind),ssnp]=="A",0,1)
  smelt$allele=alleles
  #smelt$variable=as.character(smelt$variable)
  names(smelt)=c('value','se','tvalue','variable','allele')
  #allmelt=rbind(hmelt,fmelt,smelt)

  fmelt$allele=smelt$allele
  fmelt$s_value=smelt$value
  fmelt$s_se=smelt$se
  fmelt$s_tvalue=smelt$tvalue
  #fmelt$s_value=sapply(seq(1,16), function(x) ifelse(!is.na(smelt$value[x]),smelt$value[x],0))
  #count=count+1
  #f_h_plots[[count]]=ggplot(allmelt,aes(x=variable_f,y=value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=value-se,ymax=value+se),width=.2,position=position_dodge())+facet_grid(method ~ .) + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))
  count=count+1
  s_f_h_plots[[count]]=list(values=fmelt,id=q,chrom=chr,fsnp=fsnp,hsnp=hsnp,ssnp=ssnp)
  #a<-ggplot(fmelt,aes(x=variable_f,y=f_value)) + geom_bar(stat="identity")+ geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge())+ ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))
  #b<-ggplot(fmelt,aes(x=variable_f,y=h_value)) + geom_bar(aes(fill=hapgrp),stat="identity")+ geom_errorbar(aes(ymin=h_value-h_se,ymax=h_value+h_se),width=.2,position=position_dodge())+ ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('Highest SNP effect Size of %s (%s in %s on Chr %.0f)',q,pheno,env,chr)) + theme(axis.text.x=element_text(size=8))

}

saveRDS(s_f_h_plots,'GridLMM/effect_sizes/S_F_H_all_effect_sizes.rds')
