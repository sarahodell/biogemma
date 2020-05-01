#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
e=as.character(args[[1]])


library('ggplot2')
library('data.table')
library('reshape2')
library('tibble')
library('dplyr')
library('tidyr')
library('cowplot')
# Libraries ====
library('readr')
library('ggrepel')
library('RColorBrewer')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

threshtable=fread('threshold_table.txt',data.table=F)

threshtable=threshtable[threshtable$environment==e,]
threshtable$phenotype=paste0(threshtable$phenotype,"_P")
threshtable=threshtable[threshtable$phenotype=="male_flowering_d6_P",]

gg.manhattan2 <- function(df, threshold, col, ylims){
  # format df
  df.tmp <- df %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot)) %>%
    gather(key, value, -BP,-SNP,-CHR,-BPcum,-tot)

    # Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) #%>%
    #mutate( is_annotate=ifelse(P < threshold, "yes", "no"))

  df.tmp$sig=df.tmp$value < 10^-threshold
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  ggplot(df.tmp, aes(x=BPcum, y=-log10(value))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis

    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +

    # add genome-wide sig and sugg lines
    geom_hline(yintercept = threshold,linetype="dashed") +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +

    # Add highlighted points
    geom_point(data=subset(df.tmp, sig==T), color="coral2", size=2) +

    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +

    # Custom the theme:
    theme_classic() +
    theme(
      text = element_text(size=20),
      plot.title = element_text(hjust = 0.5),
      #legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + guides(color=F)

}

snp_gwas=fread(sprintf('result_tables/600K_GWAS_%s_results.txt',e),data.table=F)
snp_gwas=snp_gwas[complete.cases(snp_gwas),]

snpthresh=threshtable[threshtable$method=="600K_SNP",]
#snpthresh=snpthresh[,c('phenotype','threshold')]
mypalette=gray.colors(5)
#mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
threshold=snpthresh[snpthresh$phenotype=="male_flowering_d6_P",]$threshold


title=sprintf("SNP GWAS %s",e)
df=snp_gwas[,c('SNP','CHR','BP','male_flowering_d6_P')]

a2<-gg.manhattan2(df,threshold,
             col=mypalette,
             ylims=c(0,10)) + labs(caption = title)

fp_gwas=fread(sprintf('result_tables/Founder_GWAS_%s_results.txt',e),data.table=F)
fpthresh=threshtable[threshtable$method=="founder_probs",]
fpthresh=fpthresh[fpthresh$phenotype=="male_flowering_d6_P",]$threshold

title=sprintf("Founder GWAS %s",e)
df=fp_gwas[,c('SNP','CHR','BP','male_flowering_d6_P')]
b2<-gg.manhattan2(df,fpthresh,
             col=mypalette,
             ylims=c(0,10)) + labs(caption = title)


hp_gwas=fread(sprintf('result_tables/Haplotype_GWAS_%s_results.txt',e),data.table=F)
hp_gwas=hp_gwas[,!names(hp_gwas) %in% "HAPGRP"]
hpthresh=threshtable[threshtable$method=="haplotype_probs",]
hpthresh=hpthresh[hpthresh$phenotype=="male_flowering_d6_P",]$threshold

title=sprintf("Haplotype GWAS %s",e)
df=hp_gwas[,c('SNP','CHR','BP','male_flowering_d6_P')]

c2<-gg.manhattan2(df,hpthresh,
             col=mypalette,
             ylims=c(0,10))+ labs(caption = title)

             prow <- plot_grid(
               a2 + theme(legend.position="none"),
               b2 + theme(legend.position="none",strip.text.x=element_blank()),
               c2 + theme(legend.position="none",axis.ticks=element_blank()),

             #  align = 'vh',
               labels = c("A", "B", "C"),
               hjust = -1,
               nrow = 3,
               ncol=1
             )

png(sprintf('result_tables/Methods_Fig3_%s_x_male_flowering_d6.png',e),width=2000,height=1500)
print(plot_grid(prow))
dev.off()
