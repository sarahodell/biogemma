#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('tidyverse')

reps=5000

ft_counts=fread('selection/haplotype_ft_gene_counts.txt',data.table=F)
ft_sum=sum(ft_counts$within_gene_count)
ft_overlap=sum(ft_counts$overlapping_gene_count)
rep_counts=c()
for(i in 1:reps){
  rep=fread(sprintf('selection/reps/haplotype_random_gene_counts_rep%.0f.txt',i),data.table=F)
  rep_counts=rbind(rep_counts,c(i,sum(rep$within_gene_count),sum(rep$overlapping_gene_count)))
}

rep_counts=as.data.frame(rep_counts)
names(rep_counts)=c('rep','gene_counts','gene_overlap')

sdp=round(sd(rep_counts$gene_counts,2))
lower=mean(rep_counts$gene_counts)-(2*sdp)
upper=mean(rep_counts$gene_counts)+(2*sdp)

dens <- density(rep_counts$gene_counts)
data <- tibble(x = dens$x, y = dens$y) %>%
  mutate(variable = case_when(
    (x >= lower & x <= upper) ~ "On",
    (x < lower | x > upper) ~ "Off",
    TRUE ~ NA_character_))


png('selection/ft_within_enrichment_density.png',width=800,height=600)
print(ggplot(data,aes(x,y)) + geom_line() + geom_area(data=filter(data,variable=="On"),fill="grey") + geom_vline(xintercept=ft_sum,color='red') + xlab("Gene Counts") + ylab("Density") + ggtitle("FT Gene Counts Within Over/Under-Represented Haplotype Blocks Over Randomized Gene Count Distribution"))
dev.off()

png('selection/ft_within_enrichment_hist.png',width=800,height=600)
print(ggplot(rep_counts,aes(x=gene_counts)) + geom_histogram(binwidth=5) + geom_vline(xintercept=ft_sum,color='red') + xlab("Gene Counts") + ylab("Density") + ggtitle("FT Gene Counts Within Over/Under-Represented Haplotype Blocks Over Randomized Gene Count Distribution"))
dev.off()

sdp=round(sd(rep_counts$gene_overlap,2))
lower=mean(rep_counts$gene_overlap)-(2*sdp)
upper=mean(rep_counts$gene_overlap)+(2*sdp)

dens <- density(rep_counts$gene_overlap)
data <- tibble(x = dens$x, y = dens$y) %>%
  mutate(variable = case_when(
    (x >= lower & x <= upper) ~ "On",
    (x < lower | x > upper) ~ "Off",
    TRUE ~ NA_character_))

png('selection/ft_overlap_enrichment_density.png',width=800,height=600)
print(ggplot(data,aes(x,y)) + geom_line() + geom_area(data=filter(data,variable=="On"),fill="grey") + geom_vline(xintercept=ft_overlap,color='red') + xlab("Gene Counts") + ylab("Density") + ggtitle("FT Gene Counts Overlapping Over/Under-Represented Haplotype Blocks Over Randomized Gene Count Distribution"))
dev.off()

png('selection/ft_overlap_enrichment_hist.png',width=800,height=600)
print(ggplot(rep_counts,aes(x=gene_overlap)) + geom_histogram(binwidth=5) + geom_vline(xintercept=ft_overlap,color='red') + xlab("Gene Counts") + ylab("Density") + ggtitle("FT Gene Counts Within Over/Under-Represented Haplotype Blocks Over Randomized Gene Count Distribution"))
dev.off()
