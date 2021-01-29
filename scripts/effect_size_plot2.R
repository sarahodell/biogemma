library('data.table')
library('ggplot2')
library('reshape2')
library('tidyverse')
library('cowplot')

s_f_h_data=readRDS('GridLMM/effect_sizes/S_F_H_all_effect_sizes.rds')
l=length(s_f_h_data)
s_f_h_plots=list()
count=0
for(n in seq(1,l)){
  #print(n)
  data=as.data.frame(s_f_h_data[[n]]$values,stringsAsFactors=F)
  na=which(is.na(data$f_value))
  if(length(na)!=0){
    data[na,]$f_value=0
    data[na,]$f_se=NA
  }
  name=s_f_h_data[[n]]$id
  #print(name)
  #ylim=max(data$f_value + f_se)
  data$variable_f=factor(data$variable,levels=data[order(data$hapgrp),]$variable_f)
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[data$variable,]
  a<-ggplot(data,aes(x=variable,y=f_value)) +
   geom_bar(aes(fill=variable),stat="identity") +
    scale_fill_manual(values=colorcodes[data$variable,]$hex_color,labels=data$variable) +
     geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge()) +
      facet_grid(.~allele) + ylab("Effect Size") + xlab("Founder") +
       ggtitle(sprintf('Founder BLUP Effect Sizes of %s',name)) +
        theme(axis.text.x=element_text(size=8))
  legend <- get_legend(
  # create some space to the left of the legend
    a + theme(legend.box.margin = margin(0, 0, 0, 20))
  )



  data0=data[data$allele==0,]
  a0<-ggplot(data0,aes(x=variable,y=f_value)) +
   geom_bar(aes(fill=variable),stat="identity") +
    scale_fill_manual(values=colorcodes[data0$variable,]$hex_color,labels=data0$variable) +
     geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge()) +
       ylab("Founder Effect Size (ggd)") + xlab("Founder")  + theme(axis.text.x=element_text(size=4),axis.title.x=element_blank(),plot.subtitle=element_text(size=8)) +
       labs(subtitle="Allele 0 (Reference) of 600K SNP") + guides(fill=F) #+  theme_classic()
  #a0
  es=unique(data1$s_value)
  se=unique(data1$s_se)
  data1=data[data$allele==1,]
  a1<-ggplot(data1,aes(x=variable,y=f_value)) +
   geom_bar(aes(fill=variable),stat="identity") +
    scale_fill_manual(values=colorcodes[data1$variable,]$hex_color,labels=data1$variable) +
     geom_errorbar(aes(ymin=f_value-f_se,ymax=f_value+f_se),width=.2,position=position_dodge()) + ylab("Founder Effect Size") +
      xlab("Founder") + theme(axis.text.x=element_text(size=4),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),plot.subtitle=element_text(size=8)) +
       labs(subtitle=sprintf("Allele 1 (Alternate) of 600K SNP (Mean Effect Size = + %.2f, SE = %.2f)",es,se)) + guides(fill=F)

  hdata = data %>% group_by(hapgrp) %>% summarize(h_value=mean(h_value),h_se=mean(h_se))
  hdata=as.data.frame(hdata,stringsAsFactors=F )
  hdata$hapgrp=as.numeric(hdata$hapgrp)
  h=max(hdata$hapgrp)
  founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
  founder_name=c(rep(founders,h))
  haplos=c()
  for(i in seq(1,h)){haplos=c(haplos,rep(i,16))}
  founder_prop=c()
  for(j in seq(1,h)){
    sub=data[data$hapgrp==j,]
    f=unique(sub$variable)
    no=length(f)
    for(k in founders){
      if(k %in% f){
        founder_prop=c(founder_prop,1/no)
      }
      else{
        founder_prop=c(founder_prop,0)
      }
    }
  }
  fprop=data.frame(founder=founder_name,hapgrp=haplos,founder_prop=founder_prop)

  h_value=c()
  h_se=c()
  for(o in seq(1,h)){
    v=hdata[hdata$hapgrp==o,]$h_value
    se=hdata[hdata$hapgrp==o,]$h_se

    h_value=c(h_value,rep(v,16))
    h_se=c(h_se,rep(se,16))
  }
  fprop$h_value=h_value
  fprop$h_se=h_se
  fprop$h_value_adj=fprop$founder_prop * fprop$h_value
  b<-ggplot(fprop,aes(x=factor(hapgrp,levels=seq(1,h)),y=h_value_adj,fill=founder)) +
   geom_bar(position="stack",stat="identity")+ scale_fill_manual(values=colorcodes[data$variable,]$hex_color,labels=data$variable)+
    ylab("Haplotype Effect Size") +
    geom_errorbar(aes(ymin=h_value-h_se,ymax=h_value+h_se),width=.02,position=position_dodge()) +
    xlab("Haplotype Groups") + ggtitle(sprintf('Haplotype BLUP Effect Sizes of %s',name)) +
     theme(axis.text.x=element_text(size=6),plot.title=element_text(size=10)) + guides(fill=F)

  p1=plot_grid(a0,a1,ncol=2,rel_widths=c(9,6)) + xlab("Founders") + labs(title=sprintf("Founder BLUP Effect Sizes of %s",name))
  p2=plot_grid(p1,b,nrow=2)
  p3=plot_grid(p2,legend,ncol=2,rel_widths = c(4,.4))
  p4=p3 + labs(title=sprintf("BLUP Effect Sizes of %s",name))
  count=count+1
  s_f_h_plots[[count]]=p4
}

pdf('GridLMM/effect_sizes/S_F_H_QTL_Effect_Sizes.pdf',onefile=TRUE)
for(i in 1:length(s_f_h_plots)){
  print(s_f_h_plots[[i]])
}
dev.off()
