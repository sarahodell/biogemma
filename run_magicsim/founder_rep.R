#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
rep=as.numeric(args[[2]])
library('data.table')
library('ggplot2')

total_tests=73248


#fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
#fprobs=fprobs[[1]]
#size=dim(fprobs[[1]])[1]
null=rep(1/16,16)

#p_multinom=sapply(seq(1,dim(fsums)[2]), function(x) dmultinom(fsums[,x],prob=null))


#df=data.table(marker=dimnames(fprobs[[1]])[[2]],p_multinom=p_multinom,pos=pmap[m,]$pos)

#fwrite(df,sprintf('selection/bg%s_multinom_results.txt',chr),row.names=F,quote=F,sep='\t')


#png(sprintf('selection/founder_probs/bg%s_multinom_scan.png',chr),width=800, height=600)
#print(ggplot(df,aes(x=pos/1e6,y=-log10(p_multinom))) + geom_point() + geom_hline(yintercept=bonf) + xlab("Position (Mb)") + ylab("-log10(P-value)") + ggtitle("Multinomial Test for Over/Under Representation of Founder Alleles"))
#dev.off()


hap_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


#fprobs=readRDS(sprintf('qtl2_files/MAGIC_DHSim_rep%.0f_c%s_genoprobs.rds',rep,chr))
fprobs=readRDS(sprintf('qtl2_files/filtered/bg%s_rep%.0f_filtered_genotype_probs.rds',chr,rep))

pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
gmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%s.csv',chr),data.table=F)

m=match(dimnames(fprobs[[1]])[[2]],pmap$marker)
markers=dimnames(fprobs[[1]])[[2]]
n_markers=dim(fprobs[[1]])[2]
size=dim(fprobs[[1]])[1]
fsums=data.frame(matrix(unlist(lapply(fprobs, function(x) round(colSums(x)))),nrow=length(fprobs),byrow=T),stringsAsFactors=F)
names(fsums)=markers
p_chi=sapply(seq(1,n_markers), function(x) chisq.test(x=round(fsums[,x]),p=null)$p.value)
#p_chi=t(p_chi)
p_chi=data.frame(marker=markers,p_chi=unlist(p_chi),stringsAsFactors=F)
rownames(p_chi)=seq(1,dim(p_chi)[1])
#names(p_chi)=founders
#upper_p_chi=sapply(seq(1,dim(fsums)[2]), function(x) pchisq(round(fsums[,x]),df=15,lower.tail=F))
#df2=data.table(marker=dimnames(fprobs)[[3]],p_chi=p_chi,pos=pmap[m,]$pos)
m=match(colnames(fprobs[[1]]),pmap$marker)
p_chi$pos=pmap[m,]$pos
p_chi$cM=gmap[match(markers,gmap$marker),]$pos
#p_chi_melt=melt(p_chi,"pos")

bonf=-log10(0.05 / total_tests)

fwrite(p_chi,sprintf('selection/bg%s_founder_chisq_rep%.0f_results.txt',chr,rep),row.names=F,quote=F,sep='\t')


png(sprintf('selection/bg%s_chi_squared_rep%.0f_scan.png',chr,rep),width=1600, height=800)
print(ggplot(p_chi,aes(x=pos/1e6,y=-log10(p_chi))) + geom_point() + geom_hline(yintercept=bonf,color="red") + xlab("Position (Mb)") + ylab("Chi-Squared -log10(P-Value)") + ggtitle(sprintf("Chi-Squared Test for Representation of Founder Alleles on Chromosome %s",chr)))
dev.off()



png(sprintf('selection/bg%s_chi_squared_scan_rep%.0f_cM.png',chr,rep),width=1600, height=800)
print(ggplot(p_chi,aes(x=cM,y=-log10(p_chi))) + geom_point() + geom_hline(yintercept=bonf,color="red") + xlab("Position (cM)") + ylab("Chi-Squared -log10(P-Value)") + ggtitle(sprintf("Chi-Squared Test for Representation of Founder Alleles on Chromosome %s",chr)))
dev.off()

#upper_p_chi=sapply(seq(1,dim(fsums)[2]), function(x) chisq.test(round(fsums[,x]),null,df=15,lower.tail=F))
#upper_p_chi=t(upper_p_chi)
#upper_p_chi=as.data.frame(upper_p_chi,stringsAsFactors=F)
#names(upper_p_chi)=founders
#upper_p_chi=sapply(seq(1,dim(fsums)[2]), function(x) pchisq(round(fsums[,x]),df=15,lower.tail=F))
#df2=data.table(marker=dimnames(fprobs)[[3]],p_chi=p_chi,pos=pmap[m,]$pos)
#m=match(colnames(fprobs[[1]]),pmap$marker)
#upper_p_chi$pos=pmap[m,]$pos
#upper_p_chi_melt=melt(upper_p_chi,"pos")

#png(sprintf('selection/bg%s_chi_squared_over_scan.png',chr),width=1600, height=800)
#print(ggplot(upper_p_chi_melt,aes(x=pos/1e6,y=-log10(value),group=variable)) + geom_point(aes(color=variable)) + geom_line(aes(color=variable)) + geom_hline(yintercept=-log10(0.05/size),color="red") + xlab("Position (Mb)") + ylab("Chi-Squared Upper Tail P-Value") + ggtitle(sprintf("Chi-Squared Test for Over-Representation of Founder Alleles on Chromosome %s",chr)))
#dev.off()
