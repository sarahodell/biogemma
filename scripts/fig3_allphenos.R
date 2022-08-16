#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
e=as.character(args[[1]])
thresh=as.numeric(args[[2]])


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
founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#col=c("Anthesis-Silking Interval"=mypalette[7],"Days to Silking" = mypalette[1],"Grain Yield"=mypalette[2],"Grain Moisture"=mypalette[3],"Days to Anthesis"=mypalette[4],"Thousand-Kernel Weight"=mypalette[5],"Plant Height"=mypalette[6])
mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
col=c("asi_P.FALSE"="black","female_flowering_d6_P.FALSE"="black","grain_yield_15_P.FALSE"="black",
"harvest_grain_moisture_P.FALSE"="black","male_flowering_d6_P.FALSE"="black",
 "tkw_15_P.FALSE"="black","total_plant_height_P.FALSE"="black","asi_P.TRUE"=mypalette[7],
 "female_flowering_d6_P.TRUE"=mypalette[1],"grain_yield_15_P.TRUE"=mypalette[2],
"harvest_grain_moisture_P.TRUE"=mypalette[3],"male_flowering_d6_P.TRUE"=mypalette[4],
"tkw_15_P.TRUE"=mypalette[5] ,"total_plant_height_P.TRUE"=mypalette[6])

labels=c("Anthesis-Silking Interval","Days to Silking","Grain Yield","Grain Moisture","Days to Anthesis","Thousand Kernel Weight","Plant Height")

thresholdtable=fread(sprintf('threshold_%.2f_table.txt',thresh),data.table=F)

thresholdtable=thresholdtable[thresholdtable$environment==e,]
thresholdtable$phenotype=paste0(thresholdtable$phenotype,"_P")
thresholdtable=thresholdtable[!(thresholdtable$phenotype %in% c("male_flowering_days_P","female_flowering_days_P")), ]

rownames(thresholdtable)=seq(1,nrow(thresholdtable))
#threshtable=threshtable[threshtable$phenotype==pheno,]

qtl_bounds=fread('Biogemma_QTL.csv',data.table=F)

phenotypes=unique(thresholdtable$phenotype)


gg.manhattan2 <- function(df, threshtable, col, ylims,legend_labels=NULL){
  df.tmp <- df %>%
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(df,., by=c("CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot)) %>%
    #mutate( is_annotate=ifelse(SNP %in% hlight, "yes", "no")) %>%
    gather(key, value, -BP,-SNP,-CHR,-BPcum,-tot)
    df.tmp = df.tmp %>% left_join(threshtable,.,by=c("phenotype"="key"))
    df.tmp$sig=df.tmp$value < 10^-df.tmp$threshold

    df.tmp$sig_f = interaction(df.tmp$phenotype,df.tmp$sig)
    axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    df.tmp=df.tmp[!is.na(df.tmp$value),]
    df.tmp=df.tmp[df.tmp$value <= quantile(df.tmp$value,0.3),]
    rownames(df.tmp)=seq(1,nrow(df.tmp))
    threshold=min(df.tmp$threshold)
    ggplot(df.tmp, aes(BPcum, -log10(value))) +
    geom_point(aes(color=sig_f), alpha=0.8, size=1.5) +
    scale_color_manual(values = col) +
    #scale_color_discrete(breaks=c(levels(df.tmp$sig_f)[8:14]),palette = col,name="Phenotype") +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,10)) + # expand=c(0,0)removes space between plot area and x axis

    # add plot and axis titles
    #ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = threshold, linetype="dashed") +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +

    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +

    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +

    # Custom the theme:
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      #legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + guides(color=F)

}

snp_gwas=fread(sprintf('result_tables/600K_GWAS_%s_results.txt',e),data.table=F)
snp_gwas=snp_gwas[complete.cases(snp_gwas),]

snpthresh=thresholdtable[thresholdtable$method=="600K_SNP",]
rownames(snpthresh)=seq(1,nrow(snpthresh))
#snpthresh=snpthresh[,c('phenotype','threshold')]
#mypalette=gray.colors(5)
#mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
#threshold=snpthresh[snpthresh$phenotype==pheno,]$threshold

title=sprintf("SNP GWAS %s",e)
df=snp_gwas[,c('SNP','CHR','BP',phenotypes)]


a2<-gg.manhattan2(df,snpthresh,
             col=col,
             ylims=c(0,10)) + labs(caption = title)



snp_bounds=qtl_bounds[qtl_bounds$Method=="600K_SNP" & qtl_bounds$Environment==e & qtl_bounds$Phenotype==p,]
for(p in phenotype)
  if(dim(snp_bounds)[1]!=0){
    row.names(snp_bounds)=seq(1,dim(snp_bounds)[1])
    for(i in 1:dim(snp_bounds)[1]){
      rowdf=snp_bounds[i,]
      start=rowdf$cum_left_bound_bp
      end=rowdf$alt_cum_right_bound_bp
      a2<-a2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill=col[p])
    }
}

#png('test.png')
#print(a2)
#dev.off()

fp_gwas=fread(sprintf('result_tables/Founder_GWAS_%s_results.txt',e),data.table=F)
fpthresh=threshtable[threshtable$method=="founder_probs",]
fpthresh=fpthresh[fpthresh$phenotype==pheno,]$threshold

title=sprintf("Founder GWAS %s %s",e,p)
df=fp_gwas[,c('SNP','CHR','BP',pheno)]
b2<-gg.manhattan2(df,fpthresh,
             col=col,
             ylims=c(0,12)) + labs(caption = title)


fp_bounds=qtl_bounds[qtl_bounds$Method=="Founder_probs" & qtl_bounds$Environment==e & qtl_bounds$Phenotype==p,]
if(dim(fp_bounds)[1]!=0){
  row.names(fp_bounds)=seq(1,dim(fp_bounds)[1])
  for(i in 1:dim(fp_bounds)[1]){
    rowdf=fp_bounds[i,]
    start=rowdf$cum_left_bound_bp
    end=rowdf$alt_cum_right_bound_bp
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
             col=col,
             ylims=c(0,12))+ labs(caption = title)

for
hp_bounds=qtl_bounds[qtl_bounds$Method=="Haplotype_probs" & qtl_bounds$Environment==e & qtl_bounds$Phenotype==p,]
if(dim(hp_bounds)[1]!=0){
  row.names(hp_bounds)=seq(1,dim(hp_bounds)[1])
  for(i in 1:dim(hp_bounds)[1]){
    rowdf=hp_bounds[i,]
    start=rowdf$cum_left_bound_bp
    end=rowdf$alt_cum_right_bound_bp
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


png(sprintf('result_tables/Methods_thresh%.2f_Fig3_%s_x_%s_allphenos.png',thresh,e),width=2000,height=1500)
print(plot_grid(prow))
dev.off()
