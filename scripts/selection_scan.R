#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')
library('ggplot2')
library('reshape2')
library('dplyr')

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


fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
gmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_gmap_c%s.csv',chr),data.table=F)

m=match(dimnames(fprobs[[1]])[[2]],pmap$marker)
markers=dimnames(fprobs[[1]])[[2]]
n_markers=dim(fprobs[[1]])[2]
size=dim(fprobs[[1]])[1]
#fsums=data.frame(matrix(unlist(lapply(fprobs, function(x) round(colSums(x)))),nrow=length(fprobs),byrow=T),stringsAsFactors=F)
fsums=c()
for(m in 1:n_markers){
  X = do.call(cbind,lapply(fprobs,function(x) x[,m]))
  frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.8])))
  fsums=rbind(fsums,(c(markers[m],as.numeric(unlist(frep2)))))
}
fsums=as.data.frame(fsums,stringsAsFactors=F)
fsums[,2:17]=apply(fsums[,2:17],MARGIN=2,function(x) as.numeric(x))

names(fsums)=c('marker',hap_founders)
rownames(fsums)=markers
p_chi=sapply(seq(1,n_markers), function(x) chisq.test(x=round(fsums[x,hap_founders]),p=null)$p.value)
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
sim_threshold=8.828254
fwrite(p_chi,sprintf('selection/founder_probs/bg%s_founder_chisq_results.txt',chr),row.names=F,quote=F,sep='\t')


png(sprintf('selection/founder_probs/bg%s_chi_squared_scan.png',chr),width=1600, height=800)
print(ggplot(p_chi,aes(x=pos/1e6,y=-log10(p_chi))) + geom_point() + geom_hline(yintercept=sim_threshold,color="red") + xlab("Position (Mb)") + ylab("Chi-Squared -log10(P-Value)") + ggtitle(sprintf("Chi-Squared Test for Representation of Founder Alleles on Chromosome %s",chr)))
dev.off()



png(sprintf('selection/founder_probs/bg%s_chi_squared_scan_cM.png',chr),width=1600, height=800)
print(ggplot(p_chi,aes(x=cM,y=-log10(p_chi))) + geom_point() + geom_hline(yintercept=sim_threshold,color="red") + xlab("Position (cM)") + ylab("Chi-Squared -log10(P-Value)") + ggtitle(sprintf("Chi-Squared Test for Representation of Founder Alleles on Chromosome %s",chr)))
dev.off()

all_p_chi=c()
for(chr in 1:10){
  p_chi=fread(sprintf('selection/founder_probs/bg%s_founder_chisq_results.txt',chr),data.table=F)
  p_chi$chr=chr
  all_p_chi=rbind(all_p_chi,p_chi)
}

p=ggplot(all_p_chi,aes(x=cM,y=-log10(p_chi))) + geom_point() +
 geom_hline(yintercept=sim_threshold,color="red") + facet_grid(.~chr) +
  xlab("Position (cM)") + ylab("Chi-Squared -log10(P-Value)") +
#  ggtitle("Chi-Squared Test for Representation of Founder Alleles") +
  theme_classic() +
  theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=18),
axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))
png('selection/all_chrom_founder_chi_sq.png',width=1500,height=1000)
print(p)
dev.off()

pos_end=c()
sig=all_p_chi[-log10(all_p_chi$p_chi)>=sim_threshold,]
rownames(sig)=seq(1,nrow(sig))
for(i in 1:nrow(sig)){
  line=sig[i,]
  chr=line$chr
  m=line$marker
  fdropped=readRDS(sprintf('genotypes/probabilities/geno_probs/dropped/bg%.0f_dropped_markers_genoprobs.rds',chr))
  dropped_markers=fdropped[[which(unlist(lapply(fdropped,function(x) x$marker==m)))]]$linked
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  sub=pmap[pmap$marker %in% dropped_markers,]
  if(dim(sub)[1]!=0){
    rownames(sub)=seq(1,dim(sub)[1])
    right=which.max(sub$pos)
    right_bound=sub[right,]$pos
    pos_end=c(pos_end,right_bound)
  }else{
    pos_end=c(pos_end,line$pos+1)
  }
}

sig$pos_end=pos_end
total=2105120000
sum(sig$size)/total * 100
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

colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
per_chr=c()

plot_list=list()
count=1
hap_founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
rownames(colorcodes)=colorcodes$founder
colorcodes=colorcodes[hap_founders,]
for(chr in 1:10){
  fprobs=readRDS(sprintf('genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  gmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_gmap_c%.0f.csv',chr),data.table=F)

  m=match(dimnames(fprobs[[1]])[[2]],pmap$marker)
  markers=dimnames(fprobs[[1]])[[2]]
  n_markers=dim(fprobs[[1]])[2]
  size=dim(fprobs[[1]])[1]
  #fsums=data.frame(matrix(unlist(lapply(fprobs, function(x) round(colSums(x)))),nrow=length(fprobs),byrow=T),stringsAsFactors=F)
  fsums=c()
  for(m in 1:n_markers){
    X = do.call(cbind,lapply(fprobs,function(x) x[,m]))
    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.8])))
    fsums=rbind(fsums,(c(markers[m],as.numeric(unlist(frep2)))))
  }
  fsums=as.data.frame(fsums,stringsAsFactors=F)
  fsums[,2:17]=apply(fsums[,2:17],MARGIN=2,function(x) as.numeric(x))
  names(fsums)=c('marker',hap_founders)
  rownames(fsums)=markers
  fsums=melt(fsums,'marker')
  m=match(fsums$marker,pmap$marker)
  fsums$pos=pmap[m,]$pos
  fsums$cM=gmap[match(markers,gmap$marker),]$pos
  fsums$chr=chr
  rep=fsums %>% group_by(variable) %>% summarize(mean_rep=mean(value))
  rep$chr=chr
  per_chr=rbind(per_chr,rep)
  fwrite(fsums,sprintf('run_magicsim/selection/chr%.0f_founder_coverage.txt',chr),row.names=F,quote=F,sep='\t')
  p=ggplot(fsums,aes(x=pos/1e6,y=value,group=variable,color=variable)) +
   geom_point() + geom_line() +
   scale_color_manual(values=colorcodes[hap_founders,]$hex_color,labels=hap_founders) +
    xlab("Position (Mb)") +
    ylab("Founder Coverage") + #ggtitle(sprintf('Chromosome %.0f',chr)) +
    theme_classic()
  plot_list[[count]]=p
  count=count+1
}
pdf('run_magicsim/selection/coverage.pdf')
for(i in seq(1,length(plot_list))){
  print(plot_list[[i]])
}
dev.off()

png('run_magicsim/per_chr_coverage.png')
print(ggplot(per_chr,aes(x=factor(chr,levels=seq(1,10)),y=mean_rep)) + geom_boxplot())
