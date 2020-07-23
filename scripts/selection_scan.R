#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])


library('data.table')
library('ggplot2')

fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs_010319.rds',chr))
fprobs=fprobs[[1]]
size=dim(fprobs)[1]
null=rep(1/16,16)
fsums=sapply(seq(1,dim(fprobs)[3]), function(j) sapply(seq(1,dim(fprobs)[2]),function(i) sum(fprobs[,i,j])))
p_multinom=sapply(seq(1,dim(fsums)[2]), function(x) dmultinom(fsums[,x],prob=null))

pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
m=match(dimnames(fprobs)[[3]],pmap$marker)
df=data.table(marker=dimnames(fprobs)[[3]],p_multinom=p_multinom,pos=pmap[m,]$pos)

fwrite(df,sprintf('selection/bg%s_multinom_results.txt',chr),row.names=F,quote=F,sep='\t')

bonf=-log10(0.05 / dim(fprobs)[3]) * 16
png(sprintf('selection/bg%s_multinom_scan.png',chr),width=800, height=600)
print(ggplot(df,aes(x=pos/1e6,y=-log10(p_multinom))) + geom_point() + geom_hline(yintercept=bonf) + xlab("Position (Mb)") + ylab("-log10(P-value)") + ggtitle("Multinomial Test for Over/Under Representation of Founder Alleles"))
dev.off()
