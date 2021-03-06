#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
e=as.character(args[[1]])
p=as.character(args[[2]])
thresh=as.numeric(args[[3]])


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
pheno=paste0(p,'_P')
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

threshtable=fread(sprintf('threshold_%.2f_table.txt',thresh),data.table=F)

threshtable=threshtable[threshtable$environment==e,]
phenotypes=unique(threshtable$phenotype)

threshtable$phenotype=paste0(threshtable$phenotype,"_P")
#threshtable=threshtable[threshtable$phenotype==pheno,]

qtl_bounds=fread('Biogemma_QTL.csv',data.table=F)

mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
pcol=list("asi_P"=mypalette[7],"female_flowering_d6_P" = mypalette[1],"grain_yield_15_P"=mypalette[2],"harvest_grain_moisture_P"=mypalette[3],"male_flowering_d6_P"=mypalette[4],"tkw_15_P"=mypalette[5],"total_plant_height_P"=mypalette[6])



gg.manhattan2 <- function(df, threshold, col, ylims,bounds){
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
  df.tmp=df.tmp[!is.na(df.tmp$value),]
  df.tmp=df.tmp[df.tmp$value <= quantile(df.tmp$value,0.3),]
  rownames(df.tmp)=seq(1,nrow(df.tmp))
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
greypalette=gray.colors(5)
#mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
threshold=snpthresh[snpthresh$phenotype==pheno,]$threshold

title=sprintf("SNP GWAS %s",e)
df=snp_gwas[,c('SNP','CHR','BP',pheno)]

a2<-gg.manhattan2(df,threshold,
             col=greypalette,
             ylims=c(0,16)) + labs(caption = title)
#qtl_bounds$height=1
#a3<-ggplot(qtl_bounds,x=alt_cum_right_bound_bp,y=height) + geom_segment(aes(x=cum_left_bound_bp,xend=alt_cum_right_bound_bp,y=height,yend=height),alpha=0.8,size=10,color=pcol[[pheno]]))

#phenotypes=c("male_flowering_d6","female_flowering_d6","total_plant_height","harvest_graint_moisture","grain_yield_15","tkw_15")

snp_bounds=qtl_bounds[qtl_bounds$Method=="600K_SNP" & qtl_bounds$Environment==e & qtl_bounds$Phenotype==p,]
if(dim(snp_bounds)[1]!=0){
  row.names(snp_bounds)=seq(1,dim(snp_bounds)[1])
  for(i in 1:dim(snp_bounds)[1]){
    rowdf=snp_bounds[i,]
    start=rowdf$cum_left_bound_bp
    end=rowdf$alt_cum_right_bound_bp
    p_P=paste0(p,'_P')
    a2<-a2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill="coral2")
  }
}
#ft_days=c("male_flowering_days","female_flowering_days")
#snp_bounds=snp_bounds[!snp_bounds$Phenotype %in% ft_days,]
#snp_bounds$p_P=paste0(snp_bounds$Phenotype,'_P')
#snp_bounds$heights=c(12,12,12.5,13,13,13.5,14)
#for(i in 1:dim(snp_bounds)[1]){
#  rowdf=snp_bounds[i,]
#  start=rowdf$cum_left_bound_bp
#  end=rowdf$alt_cum_right_bound_bp
#  ph=rowdf$Phenotype
#  h=rowdf$heights
#  p_P=paste0(ph,'_P')
#  a2<-a2+geom_segment(data=snp_bounds,aes(x=start,xend=end,y=h,yend),alpha=0.5,size=15,color=pcol[[p_P]])
#}
#a2<-a2+geom_segment(aes(x=snp_bounds$cum_left_bound_bp,xend=snp_bounds$alt_cum_right_bound_bp,y=snp_bounds$heights,yend=snp_bounds$heights+0.5),alpha=0.5,size=15,color=sapply(seq(1,nrow(snp_bounds)),function(x) pcol[[snp_bounds$p_P[x]]]))
#png('test_2.png')
#print(a2)
#dev.off()

fp_gwas=fread(sprintf('result_tables/Founder_GWAS_%s_results.txt',e),data.table=F)
fpthresh=threshtable[threshtable$method=="founder_probs",]
fpthresh=fpthresh[fpthresh$phenotype==pheno,]$threshold

title=sprintf("Founder GWAS %s %s",e,p)
df=fp_gwas[,c('SNP','CHR','BP',pheno)]
b2<-gg.manhattan2(df,fpthresh,
             col=greypalette,
             ylims=c(0,15)) + labs(caption = title)


fp_bounds=qtl_bounds[qtl_bounds$Method=="Founder_probs" & qtl_bounds$Environment==e & qtl_bounds$Phenotype==p,]
if(dim(fp_bounds)[1]!=0){
  row.names(fp_bounds)=seq(1,dim(fp_bounds)[1])
  for(i in 1:dim(fp_bounds)[1]){
    rowdf=fp_bounds[i,]
    start=rowdf$cum_left_bound_bp
    end=rowdf$alt_cum_right_bound_bp
    p_P=paste0(p,'_P')
    b2<-b2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill="coral2")
  }
}

hp_gwas=fread(sprintf('result_tables/Haplotype_GWAS_%s_results.txt',e),data.table=F)
hp_gwas=hp_gwas[,!names(hp_gwas) %in% "HAPGRP"]
hpthresh=threshtable[threshtable$method=="haplotype_probs",]
hpthresh=hpthresh[hpthresh$phenotype==pheno,]$threshold

title=sprintf("Haplotype GWAS %s %s",e,p)
df=hp_gwas[,c('SNP','CHR','BP',pheno)]



c2<-gg.manhattan2(df,hpthresh,
             col=greypalette,
             ylims=c(0,15))+ labs(caption = title)

hp_bounds=qtl_bounds[qtl_bounds$Method=="Haplotype_probs" & qtl_bounds$Environment==e & qtl_bounds$Phenotype==p,]
if(dim(hp_bounds)[1]!=0){
  row.names(hp_bounds)=seq(1,dim(hp_bounds)[1])
  for(i in 1:dim(hp_bounds)[1]){
    rowdf=hp_bounds[i,]
    start=rowdf$cum_left_bound_bp
    end=rowdf$alt_cum_right_bound_bp
    p_P=paste0(p,'_P')
    c2<-c2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill="coral2")
  }
}

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


png(sprintf('result_tables/Methods_thresh%.2f_Fig3_%s_x_%s.png',thresh,e,p),width=2000,height=1500)
print(plot_grid(prow))
dev.off()
