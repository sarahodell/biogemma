library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

s_f_h_data=readRDS('GridLMM/effect_sizes/All_effect_sizes.rds')
l=length(s_f_h_data)
s_f_h_plots=list()
acount=0
for(n in seq(1,l)){
  #print(n)
  data=as.data.frame(s_f_h_data[[n]]$values,stringsAsFactors=F)
  chr=s_f_h_data[[n]]$chrom
  grouptype=s_f_h_data[[n]]$label
  na=which(is.na(data$f_value))
  if(length(na)!=0){
    data[na,]$f_value=0
    data[na,]$f_se=NA
  }
  name=s_f_h_data[[n]]$id
  #print(name)
  #ylim=max(data$f_value + f_se)
  #data$variable_f=factor(data$variable,levels=data[order(data$hapgrp),]$variable_f)
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[data$founder,]
  leg<-ggplot(data,aes(x=variable_f,y=f_value,fill=variable_f)) +
   geom_bar(stat="identity") +
   scale_fill_manual(values=colorcodes[levels(data$variable_f),]$hex_color,labels=levels(data$variable_f)) +
   guides(fill=guide_legend(title="Founder")) +
    theme(legend.title=element_text(size=8),legend.text=element_text(size=6))


  legend <- get_legend(
     leg + theme(legend.box.margin = margin(0, 0, 0, 20))
  )

  a<-ggplot(data,aes(x=variable_f,y=f_value,color=variable_f)) +
   geom_point() +
   scale_color_manual(values=colorcodes[levels(data$variable_f),]$hex_color,labels=levels(data$variable_f))+
   geom_errorbar(aes(ymin=f_value-(2*f_se),ymax=f_value+(2*f_se)),width=.2,position=position_dodge()) +
   geom_text(aes(label=allele),vjust=-0.3,hjust=-0.3,color="black",size=3.5) +
   ylab("Effect Size") + xlab("Founder") +
   ggtitle(sprintf('%s Effect Sizes of %s',grouptype,name)) +
   theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5),
   axis.title.y=element_text(size=6),axis.title.x=element_text(size=8)) +
    guides(color=F)
  sdata = data %>% group_by(allele) %>% summarize(s_value=mean(s_value),s_se=mean(s_se))

  b<-ggplot(sdata,aes(x=factor(allele,levels=c(0,1)),y=s_value)) +
  geom_point() + geom_errorbar(aes(ymin=s_value-(2*s_se),ymax=s_value+(2*s_se))) +
  ylab("SNP Effect Size") +
  xlab("Allele") +
  #ggtitle(sprintf('SNP Effect Sizes of %s',name)) +
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5),
  axis.title.x=element_text(size=8),axis.title.y=element_text(size=6))

  hdata = data %>% group_by(hapgrp) %>% summarize(h_value=mean(h_value),h_se=mean(h_se))
  hdata$variable_h=factor(paste0("HAPGRP_",hdata$hapgrp),levels=paste0("HAPGRP_",hdata$hapgrp))
  c<-ggplot(hdata,aes(x=variable_h,y=h_value)) + geom_point() +
  geom_errorbar(aes(ymin=h_value-(2*h_se),ymax=h_value+(2*h_se))) +
  ylab("Haplotype Effect Size") +
  xlab("Haplotype Groups") +
  theme(axis.text.x=element_text(size=5),
  axis.text.y=element_text(size=5),
  axis.title.y=element_text(size=6),
  axis.title.x=element_text(size=8))

  c2=list()
  count=1
  h=max(data$hapgrp)
  for(i in seq(1,h)){
    sub=data[data$hapgrp==i,]
    c2[[count]]<- ggplot(sub,aes(x="",y=h_perc,fill=variable_f))+
     geom_bar(stat="identity",width=1) +
     coord_polar("y",start=0) +
     scale_fill_manual(values=colorcodes[sub$variable_f,]$hex_color,labels=levels(sub$variable_f)) +
     guides(fill=F) +
     labs(caption=paste0('HAPGRP_',i)) +
     theme(axis.title.x=element_blank(),
     axis.text.y=element_blank(),axis.title.y=element_blank(),
     axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),
     axis.text.x=element_blank(),plot.caption=element_text(size=5))

     count=count+1
  }

  pies=plot_grid(plotlist=c2,nrow=1)
  #png('GridLMM/effect_sizes/pie_test.png',width=1000,height=200)
  #print(p3)
  #dev.off()
  #c2<- ggplot(data,aes(x="",y=perc,fill=variable_f))+ geom_bar(stat="identity",width=1) + coord_polar("y",start=0)
  #b<-ggplot(fprop,aes(x=factor(hapgrp,levels=seq(1,h)),y=h_value_adj,fill=founder)) +
  # geom_bar(position="stack",stat="identity")+
  # scale_fill_manual(values=colorcodes[levels(fprop$founder),]$hex_color,labels=levels(fprop$founder))+
  #  ylab("Haplotype Effect Size") +
  #  geom_errorbar(aes(ymin=h_value-(2*h_se),ymax=h_value+(2*h_se)),width=.01,position=position_dodge()) +
  #  xlab("Haplotype Groups") + ggtitle(sprintf('Haplotype BLUP Effect Sizes of %s',name)) +
  #   theme(axis.text.x=element_text(size=6),plot.title=element_text(size=10)) + guides(fill=F)

  #b=b + geom_errorbar(aes(ymin=h_value-(2*h_se),ymax=h_value+(2*h_se)),width=.01,position=position_dodge()) +
  #p1=plot_grid(a0,a1,ncol=2,rel_widths=c(9,6)) + xlab("Founders") + labs(title=sprintf("Founder BLUP Effect Sizes of %s",name))
  #p2=plot_grid(p1,b,nrow=2)
  #

  p1=plot_grid(pies,c,ncol=1,rel_heights=c(3,12),rel_widths=c(1,1.4))
  p2=plot_grid(b,p1,nrow=1,rel_widths=c(3,h))
  p3=plot_grid(a,p2,nrow=2,rel_heights=c(8,8))
  #p3=plot_grid(p2,legend,ncol=2,rel_widths = c(4,.4))
  p4=plot_grid(p3,legend,ncol=2,rel_widths=c(4,.4))#+ labs(title=sprintf("Effect Sizes of %s by Method",name))

  acount=acount+1
  s_f_h_plots[[acount]]=p4
}

pdf('GridLMM/effect_sizes/All_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(s_f_h_plots)){
  print(s_f_h_plots[[i]])
}
dev.off()

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
name="male_flowering_d6_ALL_qDTA8"
chr=8
s_f_h_data=readRDS('GridLMM/effect_sizes/All_QTL_ES.rds')


n=which(unlist(lapply(s_f_h_data,function(x) x$qtl_id == name & x$focal=="Founder_probs")))

int=c(142,149,162)

snp=ses[[142]]$snp
geno=fread('')
#n=23
i=n[1]
ses=s_f_h_data
name=ses[[i]]$qtl_id
#chr=unique(qtl[qtl$pheno_env_id==name,]$Chromosome)
ssnp=ses[[i]]$keptsnp
smelt=ses[[i]]$SE
fmelt=ses[[i+1]]$SE
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
hmelt=ses[[i+2]]$SE
hsnp=ses[[i+2]]$keptsnp
h=ses[[i+2]]$hapgrp
allhaps=paste0('HAPGRP_',seq(1,h))
if(nrow(hmelt)<h){
  rown=rownames(hmelt)
  naf=allhaps[!allhaps %in% rownames(hmelt)]
  for(f in naf){
    line=rep(NA,5)
    #rownames(line)=f
    hmelt=rbind(hmelt,line)
  }
  rownames(hmelt)=c(rown,naf)
  hmelt=hmelt[allhaps,]
}
qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
pos=ses[[i]]$pos
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
focal=ses[[i]]$focal



#fdata=as.data.frame(s_f_h_data[[n]]$values,stringsAsFactors=F)
#chr=s_f_h_data[[n]]$chrom
#grouptype=s_f_h_data[[n]]$label
#na=which(is.na(data$f_value))
#if(length(na)!=0){
#  data[na,]$f_value=0
#  data[na,]$f_se=NA
#}
#name=s_f_h_data[[n]]$id
#print(name)
#ylim=max(data$f_value + f_se)
#data$variable_f=factor(data$variable,levels=data[order(data$hapgrp),]$variable_f)
colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
founderlabels=c("B73","OH43","A632","CO255","F252","A654","DK63","F2","EP1","D105","W117","F492","C103","B96","ND245","VA85")
hap_means=fmelt %>% group_by(hapgrp) %>% summarize(mean=mean(f_value))
fmelt$hap_value=hap_means[fmelt$hapgrp,]$mean
fmelt=fmelt[order(fmelt$allele,fmelt$hap_value),]
fmelt$hap_allele=paste0(fmelt$allele,'_',fmelt$hapgrp)
fmelt$variable_f=factor(fmelt$founder,levels=fmelt$founder)
colorcodes=colorcodes[fmelt$founder,]
leg<-ggplot(fmelt,aes(x=variable_f,y=f_value,fill=variable_f)) +
 geom_bar(stat="identity") +
 scale_fill_manual(values=colorcodes[levels(fmelt$variable_f),]$hex_color,labels=founder_lables) +
 guides(fill=guide_legend(title="Founder")) +
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6))


legend <- get_legend(
   leg + theme(legend.box.margin = margin(0, 0, 0, 20))
)
#fmelt$variable_f = factor(fmelt$founder,levels=c(fmelt$founder))
fmelt$has_mite = c(F,F,T,T,T,T,T,T,T,T,T,T,T,F,T,F)
label.df <- data.frame(Group = c("1_0","11_1","14_1"),
                       Value = c(25, 25,25))

                       ggplot(df, aes(x = X, y = Y)) +
                          +
                         geom_point(size = 2.5) +
                         theme_classic()

rect_left <- c(4,52,100,148)
rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 24,
  ymin = 0,
  ymax = 5)

df <- data.frame(group = c("0_8","0_5","0_5","0_3","0_12","0_1","0_1","1_13","1_14"),
                 variable_f = factor(c("EP1_inra","A654_inra","DK63","CO255_inra","F492","B73_inra","OH43_inra","ND245","VA85"),levels=levels(fmelt$variable_f)),
                 Y = c(1:9))
df$xmin=as.numeric(df$variable_f)-0.5
df$xmax=as.numeric(df$variable_f)+0.5
a=ggplot() +
geom_point(data=fmelt,aes(x=variable_f,y=f_value,color=allele)) +
  geom_rect(data=df,aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),fill="grey60",alpha=0.5)

df2=data.frame(group=c("0_8","0_9","0_5","0_6","0_3","0_10","0_12","0_2","0_1","0_4","1_13","1_7","1_14","1_11"),
  variable_f=factor(c("EP1_inra","D105_inra","A654_inra","FV2_inra","CO255_inra","W117_inra","F492","A632_usa","B73_inra","FV252_inra","ND245","C103_inra","VA85","B96"),levels=levels(fmelt$variable_f)),
  hjust=c(0.5,0.5,-2,0.5,0.5,0.5,0.5,0.5,-5,0.5,0.5,0.5,0.5,0.5),
  vjust=25)

#p + geom_text(data = label.df, label = "***")
a<-ggplot() +
geom_rect(data=df,aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),fill="grey60",alpha=0.4) +
 geom_point(data=fmelt,aes(x=variable_f,y=f_value,color=allele,group=hap_allele),position = position_dodge(width=1),size=10) +
 geom_errorbar(data=fmelt,aes(x=variable_f,ymin=f_value-(2*f_se),ymax=f_value+(2*f_se),color=allele),position = position_dodge(width=1),size=4,width=.5) +
 geom_text(data=df2,aes(x=variable_f,y=35,hjust=hjust),label=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N"),color="black",size=10)+
 geom_hline(yintercept=0,color="black") +
 ylab("DTA Effect Size (ggd)") + xlab("Founder") +
 scale_x_discrete(labels=c("EP1","D105","A654","DK63","F2","CO255","W117","F492","A632","B73","OH43","F252","ND245","C103","VA85","B96")) +
  labs(color="Allele") +
  scale_y_continuous(breaks=c(-20,-10,0,10,20,30))+
 theme_classic() +
 theme(axis.text.y=element_text(size=24),
 axis.title.y=element_text(size=30),
 axis.title.x=element_text(size=30),
 legend.title = element_text(size = 24),
 legend.text = element_text(size = 24),
 panel.grid.major.y=element_line(color="black"),
 #panel.grid.minor.y=element_line(color="black"),
  axis.text.x=element_text(size=24,vjust=0.5)) +
 theme(axis.text.x=element_text(face=ifelse(levels(fmelt$variable_f) %in% c("EP1_inra","D105_inra","A654_inra","DK63","FV2_inra","CO255_inra","W117_inra","F492","A632_usa","FV252_inra","ND245","C103_inra"),"bold.italic","plain")))


 #guides(text=F)

#  guides(color=F)

png('GridLMM/effect_sizes/vgt1_Figure4.png',width=2000,height=1000)
print(a)
dev.off()


sdata = data %>% group_by(allele) %>% summarize(s_value=mean(s_value),s_se=mean(s_se))

b<-ggplot(sdata,aes(x=factor(allele,levels=c(0,1)),y=s_value)) +
geom_point() + geom_errorbar(aes(ymin=s_value-(2*s_se),ymax=s_value+(2*s_se))) +
ylab("SNP Effect Size") +
xlab("Allele") +
#ggtitle(sprintf('SNP Effect Sizes of %s',name)) +
theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5),
axis.title.x=element_text(size=8),axis.title.y=element_text(size=6))

hdata = data %>% group_by(hapgrp) %>% summarize(h_value=mean(h_value),h_se=mean(h_se))
hdata$variable_h=factor(paste0("HAPGRP_",hdata$hapgrp),levels=paste0("HAPGRP_",hdata$hapgrp))
c<-ggplot(hdata,aes(x=variable_h,y=h_value)) + geom_point() +
geom_errorbar(aes(ymin=h_value-(2*h_se),ymax=h_value+(2*h_se))) +
ylab("Haplotype Effect Size") +
xlab("Haplotype Groups") +
theme(axis.text.x=element_text(size=5),
axis.text.y=element_text(size=5),
axis.title.y=element_text(size=6),
axis.title.x=element_text(size=8))

c2=list()
count=1
h=max(data$hapgrp)
for(i in seq(1,h)){
  sub=data[data$hapgrp==i,]
  c2[[count]]<- ggplot(sub,aes(x="",y=h_perc,fill=variable_f))+
   geom_bar(stat="identity",width=1) +
   coord_polar("y",start=0) +
   scale_fill_manual(values=colorcodes[sub$variable_f,]$hex_color,labels=levels(sub$variable_f)) +
   guides(fill=F) +
   labs(caption=paste0('HAPGRP_',i)) +
   theme(axis.title.x=element_blank(),
   axis.text.y=element_blank(),axis.title.y=element_blank(),
   axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),
   axis.text.x=element_blank(),plot.caption=element_text(size=5))

   count=count+1
}

pies=plot_grid(plotlist=c2,nrow=1)
p1=plot_grid(pies,c,ncol=1,rel_heights=c(3,12),rel_widths=c(1,1.4))
p2=plot_grid(b,p1,nrow=1,rel_widths=c(3,h))
p3=plot_grid(a,p2,nrow=2,rel_heights=c(8,8))
#p3=plot_grid(p2,legend,ncol=2,rel_widths = c(4,.4))
p4=plot_grid(p3,legend,ncol=2,rel_widths=c(4,.4))

png('GridLMM/effect_sizes/qDTA8_all_ES.png',width=1000,height=600)
print(p4)
dev.off()
