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
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)
thresh_table=fread('GridLMM/threshold_0.05_table.txt',data.table=F)
f_nots=fread('GridLMM/effect_sizes/Founder_notSNP_highest_peaks.txt',data.table=F)
s_notf=fread('GridLMM/effect_sizes/SNP_notFounder_highest_peaks.txt',data.table=F)

f_noth=fread('GridLMM/effect_sizes/Founder_notHaplotype_highest_peaks.txt',data.table=F)
h_notf=fread('GridLMM/effect_sizes/Haplotype_notFounder_highest_peaks.txt',data.table=F)

f_nots_actual=c()
for(i in 1:nrow(f_nots)){
  env=f_nots[i,]$environment
  high_snp=qtl[qtl$pheno_env_id==f_nots[i,]$pheno_env_id & qtl$Method=="Founder_probs",]$highest_SNP
  results=fread(sprintf('GridLMM/result_tables/Founder_GWAS_%s_results.txt',env),data.table=F)
  pheno=paste0(f_nots[i,]$phenotype,'_P')
  highp=-log10(results[results$SNP==high_snp,pheno])
  cutoff=thresh_table[thresh_table$phenotype==f_nots[i,]$phenotype & thresh_table$environment==env & thresh_table$method=="founder_probs",]$threshold
  f_nots_actual=c(f_nots_actual,highp-cutoff)
}
f_nots$sig_dist=f_nots_actual
names(f_nots)=c('phenotype','environment','method','chrom','highest_SNP','pheno_env_id','insig_dist','sig_dist')
fwrite(f_nots,'GridLMM/effect_sizes/Founder_notSNP_highest_peaks.txt',quote=F,row.names=F,sep='\t')
f_nots$insig_dist=-f_nots$insig_dist
f_nots_melt=melt(f_nots[,c('pheno_env_id','sig_dist','insig_dist')],'pheno_env_id')
a<-ggplot(f_nots_melt,aes(x=pheno_env_id,y=value,fill=variable)) +
 geom_bar(stat="identity",position="identity") + xlab("QTL") +
   ylab("Distance from threshold") + ggtitle("Founder not SNP QTL") + theme(axis.text.x=element_blank())

s_notf_actual=c()
for(i in 1:nrow(s_notf)){
  env=s_notf[i,]$environment
  high_snp=qtl[qtl$pheno_env_id==s_notf[i,]$pheno_env_id & qtl$Method=="600K_SNP",]$highest_SNP
  results=fread(sprintf('GridLMM/result_tables/600K_GWAS_%s_results.txt',env),data.table=F)
  pheno=paste0(s_notf[i,]$phenotype,'_P')
  highp=-log10(results[results$SNP==high_snp,pheno])
  cutoff=thresh_table[thresh_table$phenotype==s_notf[i,]$phenotype & thresh_table$environment==env & thresh_table$method=="600K_SNP",]$threshold
  s_notf_actual=c(s_notf_actual,highp-cutoff)
}
s_notf$sig_dist=s_notf_actual
names(s_notf)=c('phenotype','environment','method','chrom','highest_SNP','pheno_env_id','insig_dist','sig_dist')
fwrite(s_notf,'GridLMM/effect_sizes/SNP_notFounder_highest_peaks.txt',quote=F,row.names=F,sep='\t')
s_notf$insig_dist=-s_notf$insig_dist
s_notf_melt=melt(s_notf[,c('pheno_env_id','sig_dist','insig_dist')],'pheno_env_id')
b<-ggplot(s_notf_melt,aes(x=pheno_env_id,y=value,fill=variable)) +
 geom_bar(stat="identity",position="identity") + xlab("QTL") +
   ylab("Distance from threshold") + ggtitle("SNP not Founder QTL") + theme(axis.text.x=element_blank())

f_noth_actual=c()
for(i in 1:nrow(f_noth)){
  env=f_noth[i,]$environment
  high_snp=qtl[qtl$pheno_env_id==f_noth[i,]$pheno_env_id & qtl$Method=="Founder_probs",]$highest_SNP
  results=fread(sprintf('GridLMM/result_tables/Founder_GWAS_%s_results.txt',env),data.table=F)
  pheno=paste0(f_noth[i,]$phenotype,'_P')
  highp=-log10(results[results$SNP==high_snp,pheno])
  cutoff=thresh_table[thresh_table$phenotype==f_noth[i,]$phenotype & thresh_table$environment==env & thresh_table$method=="founder_probs",]$threshold
  f_noth_actual=c(f_noth_actual,highp-cutoff)
}
f_noth$sig_dist=f_noth_actual
names(f_noth)=c('phenotype','environment','method','chrom','highest_SNP','pheno_env_id','insig_dist','sig_dist')
fwrite(f_noth,'GridLMM/effect_sizes/Founder_notHaplotype_highest_peaks.txt',quote=F,row.names=F,sep='\t')
f_noth$insig_dist=-f_noth$insig_dist
f_noth_melt=melt(f_noth[,c('pheno_env_id','sig_dist','insig_dist')],'pheno_env_id')
c<-ggplot(f_noth_melt,aes(x=pheno_env_id,y=value,fill=variable)) +
 geom_bar(stat="identity",position="identity") + xlab("QTL") +
   ylab("Distance from threshold") + ggtitle("Founder not Haplotype QTL") + theme(axis.text.x=element_blank())


h_notf_actual=c()
for(i in 1:nrow(h_notf)){
  env=h_notf[i,]$environment
  high_snp=qtl[qtl$pheno_env_id==h_notf[i,]$pheno_env_id & qtl$Method=="Haplotype_probs",]$highest_SNP
  results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  pheno=paste0(h_notf[i,]$phenotype,'_P')
  highp=-log10(results[results$SNP==high_snp,pheno])
  cutoff=thresh_table[thresh_table$phenotype==h_notf[i,]$phenotype & thresh_table$environment==env & thresh_table$method=="haplotype_probs",]$threshold
  h_notf_actual=c(h_notf_actual,highp-cutoff)
}
h_notf$sig_dist=h_notf_actual
names(h_notf)=c('phenotype','environment','method','chrom','highest_SNP','pheno_env_id','insig_dist','sig_dist')
fwrite(h_notf,'GridLMM/effect_sizes/Haplotype_notFounder_highest_peaks.txt',quote=F,row.names=F,sep='\t')
h_notf$insig_dist=-h_notf$insig_dist
h_notf_melt=melt(h_notf[,c('pheno_env_id','sig_dist','insig_dist')],'pheno_env_id')
d<-ggplot(h_notf_melt,aes(x=pheno_env_id,y=value,fill=variable)) +
 geom_bar(stat="identity",position="identity") + xlab("QTL") +
   ylab("Distance from threshold") + ggtitle("Haplotype not Founder QTL") + theme(axis.text.x=element_blank())

png('GridLMM/result_tables/significance_distance.png',width=1000,height=1000)
print(plot_grid(a,b,c,d,ncol=2,nrow=2))
dev.off()
