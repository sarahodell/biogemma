#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

library('ggplot2')
library('data.table')

all_chroms=c()
print("Reading in files")
for(c in 1:10){
  print(c)
  mod=fread(sprintf('gemma/output/GEMMA_wgs_GWAS_%s_%s_results_chr%.0f.assoc.txt',pheno,env,c),data.table=F)
  phy_pos=unlist(lapply(strsplit(mod$rs,"_"),function(x) as.numeric(x[[2]])))

  #phy_pos=sapply(seq(1,length(mod$rs)),function(x) as.numeric(strsplit(mod$rs,'_')[[x]][[2]]))
  gwas=data.frame(chr=c,pos=phy_pos,pvalues=mod$p_wald,stringsAsFactors=F)
  all_chroms=rbind(all_chroms,gwas)
}

all_chroms$chr=as.factor(all_chroms$chr)
all_chroms$log10p = -log10(all_chroms$pvalues)
size=dim(all_chroms)[1]

#FDR Benjamini Hochberg
#all_chroms=all_chroms[order(all_chroms$pvalues),]
#rownames(all_chroms)=seq(1,size)
#all_chroms$rank=seq(1,size)
#Q=0.05
#all_chroms$qvalues=(all_chroms$rank*Q)/size
#all_chroms$log10q=-log10(all_chroms$qvalues)
#threshold=all_chroms[min(all_chroms[all_chroms$pvalues>all_chroms$qvalues,]$rank),]$log10q

print("Created all_chroms")
thresh_table=fread('GridLMM/threshold_table.txt',data.table=F,stringsAsFactors=F)
rec=thresh_table$phenotype==pheno & thresh_table$environment==env & thresh_table$method=="600K_SNP"
cutoff=thresh_table[rec,]$threshold
#cutoff=7
print(cutoff)

all_chroms$sig = all_chroms$log10p>= cutoff
# Grab only 10% lowest p-values
all_chroms=all_chroms[all_chroms$pvalues <= quantile(all_chroms$pvalues,0.1,lower.tail=T),]
label=sprintf('5%% Permutation Significance Threshold = %.2f',round(cutoff,2))

print("Plotting")
png(sprintf('gemma/images/%s_x_%s_GEMMA_WGS.png',pheno,env),width=960,height=680)
theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(size=26),axis.title=element_text(size=14,face="bold"))
theme_update(panel.background=element_blank())

print(ggplot(all_chroms,aes(x=pos/1e6,y=log10p)) + geom_point(aes(color=sig)) + scale_color_manual(breaks=all_chroms$sig,values=c("FALSE"="black","TRUE"="red")) + facet_grid(.~chr,scales="free") + ggtitle(sprintf("%s in %s Using Imputed WGS",pheno,env)) + xlab("Chromosome (Mb)") + ylab("-log10(P-Value)") + guides(color=F) + labs(caption=label))
dev.off()
