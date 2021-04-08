#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('tidyverse')
library('reshape2')
library('lme4')

chr="8"
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra",
"FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

X_full=readRDS(sprintf('../../genotypes/probabilities/geno_probs/raw/bg%s_genoprobs.rds',chr))

X_list=readRDS(sprintf('../../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
pmap=fread('../../genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)

region_start=126.5e6
region_end=126.9e6

#region_end=150.5e6
markers=pmap[pmap$pos>=region_start & pmap$pos<=region_end,]

mite_prob=fread('../mite_probabilities.txt',data.table=F)
has_mite=mite_prob[mite_prob$final >=0.9,]$ID

mnames=dimnames(X_full[[1]])[[3]]
X_region=X_full[[1]][,,mnames %in% markers$marker]
snames=dimnames(X_full[[1]])[[1]]


mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
yes_mite=founders[which(mite==T,mite)]
no_mite=founders[which(mite==F,mite)]

late_mite=yes_mite[c(1,3,5,11)]
early_mite=yes_mite[-c(1,3,5,11)]
mite_state=c("No_MITE","Late_MITE","Early_MITE","Late_MITE",
  "No_MITE","Early_MITE","Early_MITE","Late_MITE","Early_MITE",
  "Early_MITE","Early_MITE","No_MITE","Early_MITE","Late_MITE",
  "Early_MITE","No_MITE")

infolines=c()
xo_no=c()
all_breaks=c()
for (m in 1:length(snames)){
  md <- X_region[m,,] #pull out the probs for one individual
  tm <- t(md) #transpose matrix
  mdf <- as.data.frame(tm,row.names = rownames(tm))
  mdf<-rownames_to_column(mdf,"marker")
  names(mdf)=c("marker",founders)
  mdf <- merge(mdf,pmap,by.x='marker',by.y='marker') #get the physical marker positions
  mdf <- mdf[,c(2:17,19)]
  mlong<-melt(mdf,id="pos")
  #melt the dataframe
  max_prob<-mlong %>% group_by(pos) %>% top_n(1,value)
  max_prob=max_prob[max_prob$value>=0.01,] #grab the highest probability parent for each position
  max_prob<-as.data.frame(max_prob,stringsAsFactors=F)

  #t=max_prob %>% group_by(pos) %>% count

  max_prob=max_prob[max_prob$value >= 0.50,]
  max_prob<-max_prob[order(max_prob$pos),] #sort in ascending order
  rownames(max_prob)=seq(1:nrow(max_prob))
    # Find the most probable breakpoint between donors
  breaks=c()
  index=which(max_prob$variable != dplyr::lag(max_prob$variable))
  xo_no=c(xo_no,length(index))
  sample=snames[m]
  if(length(index)!=0){
    infolines=c(infolines,snames[m])
    for(i in 1:length(index)){
      chr=8
      if(i==1){
        start=region_start/1e6
        end=max_prob[index[i]-1,'pos']/1e6
        donor=as.character(max_prob[1,'variable'])
        breaks=rbind(breaks,c(sample,chr,as.numeric(start),as.numeric(end),donor,donor))
      }
      start=max_prob[index[i],'pos']/1e6
      donor=as.character(max_prob[index[i],'variable'])
      prev=i-1
      if(i==length(index)){
        end=region_end/1e6
      }
      else{
        end=max_prob[c(index[i+1]-1),'pos']/1e6
      }
      breaks=rbind(breaks,c(sample,chr,as.numeric(start),as.numeric(end),donor,donor))
    }
  }
  else{
    donor=as.character(unique(max_prob$variable))
    breaks=rbind(breaks,c(sample,chr,region_start/1e6,region_end/1e6,donor,donor))
  }
  breaks<-as.data.frame(breaks,stringsAsFactors=F)
  names(breaks)=c('sample','chr','start','end','donor1','donor2')
  all_breaks=rbind(all_breaks,breaks)
}

all_breaks$start_int=as.numeric(all_breaks$start)
all_breaks$end_int=as.numeric(all_breaks$end)
all_breaks$donor_f=factor(all_breaks$donor1,levels=founders)
all_breaks$state=mite_state[match(all_breaks$donor1,founders)]



colorcodes=fread('../effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
colorcodes=colorcodes[founders,]

left_snp=126.07
mite_start=135.946
rap27=136.009
zcn8=126.87
zcn8_end=126.883889
right_snp=150.34

sub_breaks=all_breaks[all_breaks$start_int<=zcn8 & all_breaks$end_int>=zcn8,]
fwrite(sub_breaks,'zcn8_founder_states.txt',row.names=F,quote=F,sep='\t')


sub=all_breaks[1:20,]

p <- ggplot(all_breaks,aes(x=start_int,y=1,color=donor_f)) +
    geom_segment(aes(xend=end_int,color=donor_f,yend=1),lineend="butt",size=10) +
    geom_vline(aes(xintercept=left_snp),color="black",size=1) +
    geom_vline(xintercept=mite_start,color="red",size=1) +
    geom_vline(xintercept=zcn8,color="red",size=1) +
    geom_vline(xintercept=rap27,color="red",size=1) +
    geom_vline(xintercept=right_snp,color="black",size=1) +
    facet_grid(sample ~ .) +
    ggtitle("Chromosome Breakpoints") +
    xlab("Position (Mb)")  +
    scale_color_manual(values=colorcodes[levels(all_breaks$donor_f),]$hex_color,labels=levels(all_breaks$donor_f)) +
    theme(axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=6))
#options(warn=0)
png('recombinant_qDTA8_recombinant_lines.png',width=1000,height=8000)
print(p)
dev.off()

# Grab only lines that switch between Early_MITE, Late_MITE, and No_MITE
mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
yes_mite=founders[which(mite==T,mite)]
no_mite=founders[which(mite==F,mite)]

late_mite=yes_mite[c(1,3,5,11)]
early_mite=yes_mite[-c(1,3,5,11)]

all_breaks$mite=all_breaks$donor1 %in% yes_mite
all_breaks$late_mite=all_breaks$donor1 %in% late_mite
all_breaks$state=ifelse(all_breaks$donor1 %in% no_mite, "No_MITE",ifelse(all_breaks$donor1 %in% late_mite,"Late_MITE",'Early_MITE'))

pheno="male_flowering_d6"
phenotypes=fread('../phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

all_breaks$ft_blup=data_blup[match(all_breaks$sample,data_blup$ID),]$y
all_breaks=all_breaks[order(all_breaks$ft_blup),]
all_breaks$sample_f=factor(all_breaks$sample,levels=unique(all_breaks$sample))
rownames(all_breaks)=seq(1,nrow(all_breaks))

keep=all_breaks %>% group_by(sample) %>% summarize(n=length(unique(state)))
ks=keep[keep$n>1,]$sample
sub_breaks=all_breaks[all_breaks$sample %in% ks,]

p2 <- ggplot(sub_breaks,aes(x=start_int,y=1,color=state)) +
    geom_segment(aes(xend=end_int,color=state,yend=1),lineend="butt",size=10) +
    geom_vline(aes(xintercept=left_snp),color="black",size=1) +
    geom_vline(xintercept=mite_start,color="red",size=1) +
    geom_vline(xintercept=zcn8,color="red",size=1) +
    geom_vline(xintercept=rap27,color="red",size=1) +
    geom_vline(xintercept=right_snp,color="black",size=1) +
    facet_grid(sample_f ~ .) +
    ggtitle("Chromosome Breakpoints") +
    xlab("Position (Mb)")  +
#  scale_color_manual(values=colorcodes[levels(all_breaks$donor_f),]$hex_color,labels=levels(all_breaks$donor_f)) +
    theme(axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=6))
#options(warn=0)
png('recombinant_qDTA8_recombinant_lines_by_state.png',width=4000,height=8000)
print(p2)
dev.off()

# Grab only lines that have a break between vgt1 and vgt2
k2=c()
samp=unique(sub_breaks$sample)
for(s in samp){
  t=sub_breaks[sub_breaks$sample==s,]
  for(i in seq(1,(nrow(t)-1))){
    start=t$end_int[i]
    end=t$start_int[i+1]
    if(start>=zcn8 & end<=mite_start){
      k2=c(k2,s)
    }
  }
}
k2=unique(k2)
#has_split=sub_breaks %>% group_by(sample) %>% summarize(t=sum(end>=zcn8 & start<=mite_start))
#k2=has_split[has_split$t>=1,]$sample
sub_breaks2=sub_breaks[sub_breaks$sample %in% k2,]
p3 <- ggplot(sub_breaks2,aes(x=start_int,y=1,color=state)) +
    geom_segment(aes(xend=end_int,color=state,yend=1),lineend="butt",size=20) +
    geom_vline(aes(xintercept=left_snp),color="black",size=1) +
    geom_vline(xintercept=mite_start,color="red",size=1) +
    geom_vline(xintercept=zcn8,color="red",size=1) +
    geom_vline(xintercept=rap27,color="red",size=1) +
    geom_vline(xintercept=right_snp,color="black",size=1) +
    facet_grid(sample_f ~ .) +
    ggtitle("Chromosome Breakpoints") +
    xlab("Position (Mb)")  +
#  scale_color_manual(values=colorcodes[levels(all_breaks$donor_f),]$hex_color,labels=levels(all_breaks$donor_f)) +
    theme(axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=12))
#options(warn=0)
png('recombinant_qDTA8_recombinant_lines_by_state2.png',width=6000,height=8000)
print(p3)
dev.off()

vgt1=c()
vgt2=c()
left_b=c()
right_b=c()
samp=unique(all_breaks$sample)
samp=samp[!(samp%in%c("EB.09S.H.00135","EB.09S.H.00169","EB.09S.H.00223","EB.10H.H.00143","EB.09S.H.00179","EB.10H.H.00023"))]
for(s in samp){
  t=all_breaks[all_breaks$sample==s,]
  v1=t[t$start_int<=mite_start & t$end_int>=mite_start,]$state
  vgt1=c(vgt1,v1)
  v2=t[t$start_int<=zcn8 & t$end_int>=zcn8,]$state
  vgt2=c(vgt2,v2)
  print(s)
  #print(v1)
  #print(v2)
  r=t[t$start_int<=right_snp & t$end_int>=right_snp,]$state
  right_b=c(right_b,r)
  print(r)
  l=t[t$start_int<=left_snp & t$end_int>=left_snp,]$state
  left_b=c(left_b,l)
  print(l)
}


states=data.frame(sample=samp,vgt1=factor(vgt1,levels=c("Early_MITE","Late_MITE")),vgt2=factor(vgt2,levels=c("No_MITE","Early_MITE","Late_MITE")),right_b=right_b,left_b=left_b,ft_blup=all_breaks[match(samp,all_breaks$sample),]$ft_blup,stringsAsFactors=F)
m1=lm(ft_blup ~ vgt1,data=states)

m2=lm(ft_blup ~ vgt1*right_b,data=states)

m3=lm(ft_blup ~ vgt1*left_b,data=states)
