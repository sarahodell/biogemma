#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('dplyr')
library('reshape2')
library('cowplot')
library('tibble')

#hex_colors=c("#f42896","#84ef7c","#a8bc44","#8ed1d6","#702349",
#             "#f2875b","#28ad26","#afd3ef","#937266","#56cc59",
#             "#663dd3","#478959","#47145b","#7c2126","#ad147a",
#             "#afb735")
#color schemes
#colors=c('red','grey','lightblue','yellow','blue','darkred','darkgrey','black','orange','purple','lightgreen','pink','darkgreen','darkblue','green','brown')

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
colorcodes=colorcodes[founders,]

# Chromosome eight maps
gmap=fread('../genotypes/qtl2/startfiles/Biogemma_gmap_c10.csv',data.table=F)
pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c10.csv',data.table=F)

pr10=readRDS('../genotypes/probabilities/geno_probs/bg10_filtered_genotype_probs.rds')
eight=lapply(pr10, function(x) x[240:249,])

samples=dimnames(eight[[1]])[[1]]

# In physical distance
all_pr=c()
for(m in 1:10){
  poi=lapply(eight,function(x) x[m,])
  poi=data.frame(matrix(unlist(poi),nrow=length(poi),byrow=T))
  tm=as.data.frame(t(poi),stringsAsFactors = F)
  rownames(tm)=dimnames(eight[[1]])[[2]]
  tm<-rownames_to_column(tm,"marker")
  names(tm)=c("marker",founders)
  tm <- merge(tm,pmap,by.x='marker',by.y='marker')
  tm <- tm[,c(2:17,19)]
  mlong<-melt(tm,id="pos")
  mlong$sample=samples[m]
  all_pr=rbind(all_pr,mlong)
}

all_pr$founder_f=factor(all_pr$variable,levels=founders)
use=c('EB.10H.H.00053','EB.10H.H.00052')
sub_pr=all_pr[all_pr$sample %in% use,]
e=ggplot(data=sub_pr,aes(x=pos/1e6, y=value,color=founder_f)) + facet_grid(sample~.) +
scale_color_manual(values=colorcodes$hex_color) + scale_fill_manual(values=colorcodes$hex_color) +
geom_ribbon(aes(ymin=0,ymax=value,fill=founder_f),alpha=0.5) + geom_line() +
xlab("Position (Mb)") + ylab("Probability") + scale_y_continuous(breaks=c(0.0,0.5,1.0)) +
theme_classic() + theme(strip.text.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
labs(fill="Founder") + guides(color=F)

legend <- get_legend(
  # create some space to the left of the legend
  e + theme(legend.box.margin = margin(0, 0, 0, 12))
)
e=e+guides(color=F,fill=F)

pairwise_ibd=fread('refinedibd/600K/Biogemma_600K_Founders_RefinedIBD_chr10.ibd',data.table=F)
names(pairwise_ibd)=c('ID_1','V2','ID_2','V4','CHR','start','end','V8','V9')
poi=c('F492','A654_inra')
pairwise_ibd=pairwise_ibd[pairwise_ibd$ID_1 %in% poi & pairwise_ibd$ID_2 %in% poi,]

f<-ggplot(pairwise_ibd, aes(start/1e6, 1)) +
geom_segment(aes(xend = end/1e6,yend=1),lineend="butt",size=50) + labs(caption = "Pairwise IBD Between A653 and F492 on Chromosome 10") + xlab("Position (Mb)") + ylab("Identity-By-Descent State") + xlim(c(0,150)) + ylim(c(1,1)) + theme_classic() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
ef=plot_grid(e,f,nrow=2)
ef = plot_grid(ef,legend,ncol=2,rel_widths = c(6,1.2))

all_ibd=c()
all_haps=list()
count=1
for(c in 1:10){
  ibd=fread(sprintf('refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',c),data.table=F)
  #all_haps=rbind(all_haps,ibd)

  p = ggplot(ibd,aes(ymin=0)) + geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymax = n_haps)) + labs(subtitle=c) + theme_classic()
  #if(c==5){
  #  p = p + xlab("Position (Mb)") + theme(axis.title.y=element_blank())
  #}

  if(c==1){
    p = p + ylab("Number of Haplotypes") + theme(axis.title.x=element_blank())
  }
  else{
    p = p + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  }
  all_haps[[count]]=p
  all_ibd=rbind(all_ibd,ibd)
  count=count+1
}

allp=plot_grid(plotlist=all_haps,ncol=10)

png('no_haps_per_chrom.png',width=2000,height=800)
print(allp)
dev.off()

all_ibd$size=all_ibd$end-all_ibd$start
h=ggplot(all_ibd,aes(x=as.factor(n_haps),y=log10(size))) + geom_jitter(alpha=0.6) + geom_boxplot(alpha=0.75,fill='grey') + theme_classic() + xlab('Number of Unique Haplotypes') + ylab('Haplotype Block Size log10(bp)')

png('hap_size_dist.png',width=800,height=800)
print(h)
dev.off()

p1=plot_grid(ef,h,ncol=2,rel_widths=c(4,3))
p2=plot_grid(allp,p1,nrow=2,rel_heights=c(6,4),labels = c("A", "B", "C"))

png('figure2.png',width=2000,height=1000)
print(p2)
dev.off()
