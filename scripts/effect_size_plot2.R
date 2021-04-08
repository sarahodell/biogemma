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

png('GridLMM/effect_sizes/test.png',width=1000,height=600)
print(p4)
dev.off()
