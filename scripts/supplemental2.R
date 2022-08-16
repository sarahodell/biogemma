#!/usr/bin/env Rscript

# Libraries ====
library('ggplot2')
library('data.table')
library('reshape2')
library('tibble')
library('tidyr')
library('dplyr')
library('cowplot')
library('readr')
library('ggrepel')
library('RColorBrewer')

plotlist=list()

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
thresh=0.05

#phenotypes=unique(threshtable$phenotype)
phenotypes=c("male_flowering_d6","female_flowering_d6",
'grain_yield_15','tkw_15','total_plant_height',"harvest_grain_moisture")
environments=c('ALL','BLOIS_2014_OPT','BLOIS_2017_OPT',
'GRANEROS_2015_OPT','NERAC_2016_WD','STPAUL_2017_WD','SZEGED_2017_OPT')

#threshtable=threshtable[threshtable$phenotype==pheno,]

bounds=fread('Biogemma_QTL_all_thresholds.csv',data.table=F)
bounds=bounds[bounds$`10P_Sig`!="T",]
pheno_envs=unique(bounds$pheno_env)
pheno_envs=pheno_envs[pheno_envs != "male_flowering_d6_ALL"]

mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
col=c("asi_P.FALSE"="black","female_flowering_d6_P.FALSE"="black","grain_yield_15_P.FALSE"="black",
"harvest_grain_moisture_P.FALSE"="black","male_flowering_d6_P.FALSE"="black",
 "tkw_15_P.FALSE"="black","total_plant_height_P.FALSE"="black","asi_P.TRUE"=mypalette[1],
 "female_flowering_d6_P.TRUE"=mypalette[2],"grain_yield_15_P.TRUE"=mypalette[7],
"harvest_grain_moisture_P.TRUE"=mypalette[6],"male_flowering_d6_P.TRUE"=mypalette[5],
"tkw_15_P.TRUE"=mypalette[3] ,"total_plant_height_P.TRUE"=mypalette[4])

gg.manhattan2 <- function(df, threshold, col, ylims,legend_labels=NULL){
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
    #df.tmp = df.tmp %>% left_join(threshtable,.,by=c("phenotype"="key"))
    #threshold=thresht
    df.tmp$sig=-log10(df.tmp$value) >= threshold

    df.tmp$sig_f = interaction(df.tmp$key,df.tmp$sig)
    axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    df.tmp=df.tmp[!is.na(df.tmp$value),]
    #df.tmp=df.tmp[df.tmp$value <= quantile(df.tmp$value,0.3),]
    rownames(df.tmp)=seq(1,nrow(df.tmp))
    theme_set(theme_classic())
    theme_update(text=element_text(family="Helvetica"))
    theme_update(plot.caption = element_text(hjust = 0))
    theme_update(plot.title = element_text(size=10),axis.title=element_text(size=10))
    theme_update(panel.background=element_blank())
    theme_update(plot.caption=element_text(size=10))
    theme_update(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))

    ggplot(df.tmp, aes(BPcum, -log10(value))) +
    geom_point(aes(color=sig_f), alpha=0.8, size=1.0) +
    scale_color_manual(values = col) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,12)) +
    labs(x = "Chromosome") +
    geom_hline(yintercept = threshold, linetype="dashed") +
    theme(
      plot.subtitle = element_text(size=10),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + guides(color=F)

}

count=1

for(p in phenotypes){
  for(e in environments){
    pe=paste0(p,'_',e)
    if(pe %in% pheno_envs){
      threshtable=fread(sprintf('threshold_%.2f_table.txt',thresh),data.table=F)
      #threshtable$phenotype=paste0(threshtable$phenotype,"_P")
      threshtable=threshtable[threshtable$environment==e,]
      pheno=paste0(p,'_P')
      snp_gwas=fread(sprintf('result_tables/600K_GWAS_%s_results.txt',e),data.table=F)
      snp_gwas=snp_gwas[complete.cases(snp_gwas),]

      snpthresh=threshtable[threshtable$method=="600K_SNP",]
      greypalette=gray.colors(5)
      snpthresh=snpthresh[snpthresh$phenotype==p,]$threshold

      title=sprintf("GWAS SNP %s %s",e,p)
      df=snp_gwas[,c('SNP','CHR','BP',pheno)]
      #cutoff=quantile(df[,pheno],0.10)
      #df=df[df[,pheno]<=cutoff,]
      df$sig=-log10(df[,pheno])>=snpthresh
      df1=df[df$sig==T,]
      df2=df[df$sig==F,]
      rdraw=sample(rownames(df2),4000)
      df2=df2[rdraw,]
      if(dim(df1)[1]!=0){
        df=rbind(df1,df2)
      }else{
        df=df2
      }
      a2<-gg.manhattan2(df,snpthresh,
               col=col,
               ylims=c(0,16)) + labs(caption = title)

      snp_bounds=bounds[bounds$Method=="600K_SNP" & bounds$Environment==e & bounds$Phenotype==p & bounds$`10P_Sig`!=T,]
      if(dim(snp_bounds)[1]!=0){
        row.names(snp_bounds)=seq(1,dim(snp_bounds)[1])
        for(i in 1:dim(snp_bounds)[1]){
          rowdf=snp_bounds[i,]
          start=rowdf$cum_left_bound_bp
          end=rowdf$cum_right_bound_bp
          color=col[paste0(pheno,'.TRUE')]
          a2<-a2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill=color)
          }
      }
      fp_gwas=fread(sprintf('result_tables/Founder_GWAS_%s_results.txt',e),data.table=F)
      fpthresh=threshtable[threshtable$method=="founder_probs",]
      fpthresh=fpthresh[fpthresh$phenotype==p,]$threshold

      title=sprintf("QTL F %s %s",e,p)
      df=fp_gwas[,c('SNP','CHR','BP',pheno)]
      b2<-gg.manhattan2(df,fpthresh,
               col=col,
               ylims=c(0,16)) + labs(caption = title)

      fp_bounds=bounds[bounds$Method=="Founder_probs" & bounds$Environment==e & bounds$Phenotype==p & bounds$`10P_Sig`!=T,]
      if(dim(fp_bounds)[1]!=0){
        row.names(fp_bounds)=seq(1,dim(fp_bounds)[1])
        for(i in 1:dim(fp_bounds)[1]){
          rowdf=fp_bounds[i,]
          start=rowdf$cum_left_bound_bp
          end=rowdf$cum_right_bound_bp
          color=col[paste0(pheno,'.TRUE')]
          b2<-b2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill=color)
        }
      }

      hp_gwas=fread(sprintf('result_tables/Haplotype_GWAS_%s_results.txt',e),data.table=F)
      hp_gwas=hp_gwas[,!names(hp_gwas) %in% "HAPGRP"]
      hpthresh=threshtable[threshtable$method=="haplotype_probs",]
      hpthresh=hpthresh[hpthresh$phenotype==p,]$threshold

      title=sprintf("QTL H %s %s",e,p)
      df=hp_gwas[,c('SNP','CHR','BP',pheno)]
      df$sig=-log10(df[,pheno])>=hpthresh
      df1=df[df$sig==T,]
      df2=df[df$sig==F,]
      rdraw=sample(rownames(df2),4000)
      df2=df2[rdraw,]
      if(dim(df1)[1]!=0){
        df=rbind(df1,df2)
      }else{
        df=df2
      }


      c2<-gg.manhattan2(df,hpthresh,
               col=col,
               ylims=c(0,16))+ labs(caption = title)

      hp_bounds=bounds[bounds$Method=="Haplotype_probs" & bounds$Environment==e & bounds$Phenotype==p & bounds$`10P_Sig`!=T,]
      if(dim(hp_bounds)[1]!=0){
        row.names(hp_bounds)=seq(1,dim(hp_bounds)[1])
        for(i in 1:dim(hp_bounds)[1]){
          rowdf=hp_bounds[i,]
          start=rowdf$cum_left_bound_bp
          end=rowdf$cum_right_bound_bp
          color=col[paste0(pheno,'.TRUE')]
          c2<-c2+geom_ribbon(aes_string(xmin=start,xmax=end),alpha=0.2,fill=color)
        }
      }

      prow2 <- plot_grid(
          a2 + theme(legend.position="none"),
          b2 + theme(legend.position="none",strip.text.x=element_blank()),
          c2 + theme(legend.position="none",axis.ticks=element_blank()),
          #  align = 'vh',
          labels = c("A", "B", "C"),
          label_size = 12,
          nrow = 3,
          ncol=1
      )
      plotlist[[count]]=prow2
      count=count+1
    }
  }
}

pdf('result_tables/Supplemental_File2.pdf')
for(i in 1:length(plotlist)){
  print(plotlist[[i]])
}
dev.off()
