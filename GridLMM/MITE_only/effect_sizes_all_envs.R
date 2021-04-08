#!/usr/bin/env Rscript
#args=commandArgs(trailingOnly=T)
#c=as.character(args[[1]])
#pheno=as.character(args[[2]])

library('data.table')
library('ggplot2')
library('dplyr')
library('reshape2')
library('lme4')
library('lme4qtl')
#library('emmeans')

c="8"
pheno="male_flowering_d6"
environments=c("ALL","STPAUL_2017_WD")
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)

colorcodes=fread('../effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
colorcodes=colorcodes[founders,]
#ibd=fread(sprintf('../../ibd_segments/bg%s_ibd_blocks.txt',c),data.table=F)
pmap=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)
bounds=fread('../Biogemma_QTL.csv',data.table=F)
threshold_table=fread('../threshold_0.05_table.txt',data.table=F)
ses=list()
count=1
#start=100*1e6
#mite_end=135.946644

#low_rep=readRDS(sprintf('../../genotypes/probabilities/geno_probs/dropped/bg%s_low_rep_markers_MITE_only.rds',c))
mite_prob=fread('../mite_probabilities.txt',data.table=F)
has_mite=mite_prob[mite_prob$final>=0.9,]$ID
for(env in environments){
  model=readRDS(sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs_MITE_only.rds',c,pheno,env))
  model_merge=merge(model,pmap,by.x="X_ID",by.y="marker")
  model_merge=model_merge[order(model_merge$pos),]
  row.names(model_merge)=seq(1,dim(model_merge)[1])
  betas=sprintf('beta.%.0f',seq(1,16))
  # Discount effect sizes of low representation founder markers
  #for(f in 1:16){
  #  low_markers=low_rep[[f]]
  #  if(!is.null(low_markers)){
  #    for(m in low_markers){
  #      model_merge[model_merge$X_ID==m,paste0('beta.',f)]=NA
  #    }
  #  }
  #}
  cutoff=threshold_table[threshold_table$phenotype==pheno & threshold_table$method=="founder_probs" & threshold_table$environment==env,]$threshold
  model_merge=model_merge[-log10(model_merge$p_value_ML)>=cutoff,]
  row.names(model_merge)=seq(1,dim(model_merge)[1])

  markers=model_merge$X_ID
  founder_probs = readRDS(sprintf('../../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',c))
  K = fread(sprintf('../K_matrices/K_matrix_chr%s.txt',c),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('../phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))

  for(m in markers){
    mpos=pmap[pmap$marker==m,]$pos
    X = do.call(cbind,lapply(founder_probs,function(x) x[has_mite,m]))
    colnames(X) = founders
    rownames(X) = has_mite
    if(env =="ALL"){
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      data_blup=data_blup[data_blup$ID %in% has_mite,]
      X = X[rownames(X) %in% data_blup$ID,]
      subK=K[rownames(X),rownames(X)]
      data_blup=data_blup[rownames(X),]

      m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      rownames(se4)=founders
      names(se4)=c('value','se','tvalue')
      se4$founder=rownames(se4)
      se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('%s_BLUP_MITE_only_founder_effect_sizes_lme4qtl.png',m),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
      theme(axis.text.x=element_text(size=10)) +
      xlab("Founder") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates Chr 8:%.0f",m,mpos)))
      dev.off()
    }
    else{
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      phenotype=phenotype[phenotype$Genotype_code %in% has_mite,]
      X = X[rownames(X) %in% phenotype$Genotype_code,]
      phenotype=phenotype[rownames(X),]
      subK=K[rownames(X),rownames(X)]
      m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=founders
      se4$founder=rownames(se4)
      se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('%s_%s_MITE_only_founder_effect_sizes_lme4qtl.png',m,env),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
       geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
       scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
      theme(axis.text.x=element_text(size=10)) +
      xlab("Founder") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates Chr 8:%.0f",m,mpos)))
      dev.off()
    }
  #sub = model_merge[,c("X_ID",betas)]
  #names(sub)=c("X_ID",founders)
  #mlong=melt(sub,'X_ID')
  #mlong$pos=model_merge[match(mlong$X_ID,model_merge$X_ID),]$pos
  #mlong$env=env
  #
  #markers=unique(mlong$X_ID)
  #for(m in markers){
  #  count=count+1
  #  eff=mlong[mlong$X_ID==m,]
  #  eff$has_mite=has_mite
  #  eff=eff[order(eff$value),]
  #  rownames(eff)=seq(1,nrow(eff))
  #  eff$founder_f=factor(eff$variable,levels=eff$variable)
  #  p<-ggplot(eff,aes(x=founder_f,y=value,color=has_mite)) + geom_bar(aes(fill=has_mite),stat="identity") + ylab("Effect Size") + xlab("Founder") + ggtitle(sprintf('%s in %s of MITE Only MAGIC Lines Chr %s : %.0f Marker %s',pheno,env,c,unique(eff$pos),m))
  #  plots[[count]]=p
  #}
  ses[[count]]=list(SE=se4,qtl_id=sprintf("male_flowering_d6_%s_qDTA8_MITE_only",env),snp=m,pos=mpos)
  count=count+1
  }
}

saveRDS(ses,'MITE_Only_Founder_prob_QTL_SEs.rds')
#mlong_all=as.data.frame(mlong_all,stringsAsFactors=F)
#names(mlong_all)=c("marker","variable","value","pos","env")

#pdf(sprintf('images/chr%s_%s_founderprobs_avg_effect_sizes_all_founders_MITE_only.pdf',c,pheno),width=15,height=15)
#for( i in 1:length(plots)){
#  print(plots[[i]])
#}
#dev.off()

#png(sprintf('images/chr%s_%s__founderprobs_avg_effect_sizes_all_founders_MITE_only.png',c,pheno,env),width=960,height=960)
#print(ggplot(mlong_all,aes(x=factor(env),y=value)) + geom_boxplot(aes(fill=env)) + facet_wrap(~variable) +
#ggtitle(sprintf("Founder Probability %s  Average vgt1 Effect Sizes by Environment (within vgt1 QTL bound for chromosome %s)",pheno,c)) + ylab("Effect Size") + geom_hline(yintercept=0,color='red') +
# xlab("Environment") +
# theme_classic() +
# theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major.y = element_line(color = "grey80"),panel.grid.minor.y = element_line(color = "grey80")))
#dev.off()

#png(sprintf('images/chr%s_%s__founderprobs_highest_SNP_effect_sizes_all_founders_MITE_only.png',c,pheno,env),width=960,height=960)
#print(ggplot(mlong_point,aes(x=factor(env),y=value)) + geom_bar(aes(fill=env),stat="identity") + facet_wrap(~variable) +
#ggtitle(sprintf("Founder Probability %s vgt1 Effect Sizes by Environment at highest SNP, chromosome %s MITE Only)",pheno,c)) + ylab("Effect Size") + geom_hline(yintercept=0,color='red') +
# xlab("Environment") +
# theme_classic() +
# theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major.y = element_line(color = "grey80"),panel.grid.minor.y = element_line(color = "grey80")))
#dev.off()
