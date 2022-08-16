#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('tidyverse')
library('cowplot')

reps=1000

ft_counts=fread('selection/interchrom_ld_ft_gene_counts.txt',data.table=F)
ft_overlap=sum(ft_counts$overlapping_gene_count)
rep_counts=c()
for(i in 1:reps){
  rep=fread(sprintf('selection/reps/interchrom_ld_random_gene_counts_rep%.0f.txt',i),data.table=F)
  rep_counts=rbind(rep_counts,c(i,sum(rep$overlapping_gene_count)))
}

rep_counts=as.data.frame(rep_counts)
names(rep_counts)=c('rep','gene_overlap')

sdp=round(sd(rep_counts$gene_overlap,2))
lower=mean(rep_counts$gene_overlap)-(2*sdp)
upper=mean(rep_counts$gene_overlap)+(2*sdp)

dens <- density(rep_counts$gene_overlap)
data <- tibble(x = dens$x, y = dens$y) %>%
  mutate(variable = case_when(
    (x >= lower & x <= upper) ~ "On",
    (x < lower | x > upper) ~ "Off",
    TRUE ~ NA_character_))

theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(plot.title = element_text(size=30),axis.title=element_text(size=24))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=30))
theme_update(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))

a<-ggdraw() + draw_image('stats/ld_decay/circos/circos-physical.png')

p=ggplot(data,aes(x,y)) + geom_line() +
 geom_area(data=filter(data,variable=="On"),fill="grey") +
  geom_vline(xintercept=ft_overlap,color='red') +
   xlab("Gene Counts") + ylab("Density") + geom_text(label=ft_overlap,x=ft_overlap-6,y=0.025,color="black",size=8)

prow=plot_grid(a,p,ncol=1,rel_heights=c(3,2),labels=c("A","B"),label_size=30)

png('selection/ft_interchrom_ld_figure.png',width=1000,height=2000)
print(prow)
dev.off()


png('selection/ft_interchrom_ld_overlap_enrichment_density.png',width=1000,height=800)
print(p)
dev.off()

png('selection/ft_interchrom_ld_overlap_enrichment_hist.png',width=800,height=600)
print(ggplot(rep_counts,aes(x=gene_overlap)) + geom_histogram(binwidth=5) + geom_vline(xintercept=ft_overlap,color='red') +
 xlab("Gene Counts") + ylab("Density") +
  ggtitle("FT Gene Counts Overlapping Regions of High Inter-chromosomal LD Over Randomized Gene Count Distribution"))
dev.off()
