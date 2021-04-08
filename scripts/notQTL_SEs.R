library('data.table')
library('lme4')
library('lme4qtl')
library('ggplot2')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)

qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
ft_days=c('male_flowering_days','female_flowering_days')
qtl_overlap=qtl_overlap[!(qtl_overlap$Phenotype %in% ft_days),]
rownames(qtl_overlap)=seq(1,nrow(qtl_overlap))

h_notsf=fread('GridLMM/result_tables/H_only_not_S_and_F_highest_peaks.txt',data.table=F)
f_notsh=fread('GridLMM/result_tables/F_only_not_S_and_H_highest_peaks.txt',data.table=F)
s_notfh=fread('GridLMM/result_tables/SNP_only_not_F_and_H_highest_peaks.txt',data.table=F)
noth=fread('GridLMM/result_tables/S_and_F_not_H_highest_peaks.txt')
nots=fread('GridLMM/result_tables/F_and_H_not_S_highest_peaks.txt')

s_f_h_ids=qtl_overlap[qtl_overlap$label=="S_F_H",]$pheno_env_id
s_only_ids=qtl_overlap[qtl_overlap$label=="S_only",]$pheno_env_id
f_only_ids=qtl_overlap[qtl_overlap$label=="F_only",]$pheno_env_id
h_only_ids=qtl_overlap[qtl_overlap$label=="H_only",]$pheno_env_id
s_and_f_ids=qtl_overlap[qtl_overlap$label=="S_and_F",]$pheno_env_id
f_and_h_ids=qtl_overlap[qtl_overlap$label=="F_and_H",]$pheno_env_id
#m1=c('Founder_probs','600K_SNP')

# S only, not F and H
ses=list()
count=1
q1=s_notfh[s_notfh$pheno_env_id %in% s_only_ids,]
rownames(q1)=seq(1,nrow(q1))
for(q in 1:nrow(q1)){
  line=q1[q,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  snp=line$highest_SNP
  name=line$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  #Founder
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[founders,]
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  if(method=="Founder_probs"){
    founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
    X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
    colnames(X) = founders
    rownames(X) = dimnames(founder_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=founders
      se4$founder=rownames(se4)
      se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_S_only_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
       geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
       scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
      theme(axis.text.x=element_text(size=10)) +
      xlab("Founder") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]

      m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      rownames(se4)=founders
      names(se4)=c('value','se','tvalue')
      se4$founder=rownames(se4)
      se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_S_only_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
      theme(axis.text.x=element_text(size=10)) +
      xlab("Founder") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  if(method=="Haplotype_probs"){
    results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
    h=results[results$SNP==snp,]$HAPGRP
    haplo_probs = readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
    X = do.call(cbind,lapply(haplo_probs,function(x) x[,snp]))
    colnames(X) = paste0("HAPGRP_",seq(1,h))
    rownames(X) = dimnames(haplo_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=colnames(X)
      se4$hapgrp=rownames(se4)
      se4$variable_f=factor(se4$hapgrp,levels=se4$hapgrp)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_S_only_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Haplotype") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }
    else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=colnames(X)
      se4$hapgrp=rownames(se4)
      se4$variable_f=factor(se4$hapgrp,levels=se4$hapgrp)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_S_only_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Haplotype") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  ses[[count]]=list(SE=se4,method=method,qtl_id=name,snp=snp,pos=pos)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/founder_ES/S_only_not_F_and_H_QTL_SEs.rds')

# F only, not S and H
ses=list()
count=1
q1=f_notsh[f_notsh$pheno_env_id %in% f_only_ids,]
rownames(q1)=seq(1,nrow(q1))
for(q in 1:nrow(q1)){
  line=q1[q,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  snp=line$highest_SNP
  name=line$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  #Founder
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[founders,]
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  if(method=="600K_SNP"){
    X = fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
    rownames(X)=X$ind
    X = X[,snp,drop=F]
    #colnames(X) = founders
    #rownames(X) = dimnames(founder_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 1 + unlist(X) + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      se4[2,]$value=se4$value[-1]+se4$value[1]
      se4$variable_f=factor(c(0,1),levels=c(0,1))
      #se4$founder=rownames(se4)
      #se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_F_only_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value)) +
      geom_point() +
       geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Allele") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]

      m4 = relmatLmer(y ~ 1 + unlist(X) + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      #rownames(se4)=founders
      names(se4)=c('value','se','tvalue')
      se4[2,]$value=se4$value[-1]+se4$value[1]
      se4$variable_f=factor(c(0,1),levels=c(0,1))

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_F_only_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Allele") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  if(method=="Haplotype_probs"){
    results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
    h=results[results$SNP==snp,]$HAPGRP
    haplo_probs = readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
    X = do.call(cbind,lapply(haplo_probs,function(x) x[,snp]))
    colnames(X) = paste0("HAPGRP_",seq(1,h))
    rownames(X) = dimnames(haplo_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=colnames(X)
      se4$hapgrp=rownames(se4)
      se4$variable_f=factor(se4$hapgrp,levels=se4$hapgrp)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_F_only_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Haplotype") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }
    else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=colnames(X)
      se4$hapgrp=rownames(se4)
      se4$variable_f=factor(se4$hapgrp,levels=se4$hapgrp)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_F_only_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Haplotype") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  ses[[count]]=list(SE=se4,method=method,qtl_id=name,snp=snp,pos=pos)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/founder_ES/F_only_not_S_and_H_QTL_SEs.rds')

# H only, not S or F
ses=list()
count=1
q1=h_notsf[h_notsf$pheno_env_id %in% h_only_ids,]
rownames(q1)=seq(1,nrow(q1))
for(q in 1:nrow(q1)){
  line=q1[q,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  snp=line$highest_SNP
  name=line$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  #Founder
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  rownames(colorcodes)=colorcodes$founder
  colorcodes=colorcodes[founders,]
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  if(method=="600K_SNP"){
    X = fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
    rownames(X)=X$ind
    X = X[,snp,drop=F]
    #colnames(X) = founders
    #rownames(X) = dimnames(founder_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 1 + unlist(X) + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      se4[2,]$value=se4$value[-1]+se4$value[1]
      se4$variable_f=factor(c(0,1),levels=c(0,1))
      #se4$founder=rownames(se4)
      #se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_H_only_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value)) +
      geom_point() +
       geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Allele") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]

      m4 = relmatLmer(y ~ 1 + unlist(X) + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      #rownames(se4)=founders
      names(se4)=c('value','se','tvalue')
      se4[2,]$value=se4$value[-1]+se4$value[1]
      se4$variable_f=factor(c(0,1),levels=c(0,1))

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_H_only_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Allele") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  if(method=="Founder_probs"){
    founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
    X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
    colnames(X) = founders
    rownames(X) = dimnames(founder_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=founders
      se4$founder=rownames(se4)
      se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_H_only_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
       geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
       scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
      theme(axis.text.x=element_text(size=10)) +
      xlab("Founder") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]

      m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      rownames(se4)=founders
      names(se4)=c('value','se','tvalue')
      se4$founder=rownames(se4)
      se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_H_only_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
      theme(axis.text.x=element_text(size=10)) +
      xlab("Founder") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  ses[[count]]=list(SE=se4,method=method,qtl_id=name,snp=snp,pos=pos)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/founder_ES/H_only_not_S_and_F_QTL_SEs.rds')

# S and F , not H
ses=list()
count=1
q1=noth[noth$pheno_env_id %in% s_and_f_ids,]
rownames(q1)=seq(1,nrow(q1))
for(q in 1:nrow(q1)){
  line=q1[q,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  snp=line$highest_SNP
  name=line$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  #Founder
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  #colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  #rownames(colorcodes)=colorcodes$founder
  #colorcodes=colorcodes[founders,]
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  if(method=="Haplotype_probs"){
    results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
    h=results[results$SNP==snp,]$HAPGRP
    haplo_probs = readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
    X = do.call(cbind,lapply(haplo_probs,function(x) x[,snp]))
    colnames(X) = paste0("HAPGRP_",seq(1,h))
    rownames(X) = dimnames(haplo_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=colnames(X)
      se4$hapgrp=rownames(se4)
      se4$variable_f=factor(se4$hapgrp,levels=se4$hapgrp)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_S_and_F_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Haplotype") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }
    else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      rownames(se4)=colnames(X)
      se4$hapgrp=rownames(se4)
      se4$variable_f=factor(se4$hapgrp,levels=se4$hapgrp)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_S_and_F_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Haplotype") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  ses[[count]]=list(SE=se4,method=method,qtl_id=name,snp=snp,pos=pos)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/founder_ES/S_and_F_not_H_QTL_SEs.rds')

# Not S, F and H
ses=list()
count=1
q1=nots[nots$pheno_env_id %in% f_and_h_ids,]
rownames(q1)=seq(1,nrow(q1))
for(q in 1:nrow(q1)){
  line=q1[q,]
  method=line$method
  pheno=line$phenotype
  env=line$environment
  chr=as.character(line$chrom)
  snp=line$highest_SNP
  name=line$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  #Founder
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  #colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
  #rownames(colorcodes)=colorcodes$founder
  #colorcodes=colorcodes[founders,]
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  if(method=="600K_SNP"){
    X = fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
    rownames(X)=X$ind
    X = X[,snp,drop=F]
    #colnames(X) = founders
    #rownames(X) = dimnames(founder_probs[[1]])[[1]]
    if(env != "ALL"){
      phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
      phenotype = phenotype[!is.na(phenotype$y),]
      phenotype$y=phenotype$y-mean(phenotype$y)
      rownames(phenotype)=phenotype$Genotype_code
      i=intersect(phenotype$Genotype_code,rownames(X))
      X = X[i,]
      phenotype=phenotype[i,]
      subK=K[i,i]
      m4 = relmatLmer(y ~ 1 + unlist(X) + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      names(se4)=c('value','se','tvalue')
      se4[2,]$value=se4$value[-1]+se4$value[1]
      se4$variable_f=factor(c(0,1),levels=c(0,1))
      #se4$founder=rownames(se4)
      #se4$variable_f=factor(se4$founder,levels=se4$founder)

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_F_and_H_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value)) +
      geom_point() +
       geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Allele") + ylab("Effect Size") +
      labs(title=sprintf("%s Effect Size Estimates",name)))
      dev.off()
    }else{
      phenotype = phenotype[!is.na(phenotype$y),]
      m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
      data_blup = as.data.frame(ranef(m0)$Genotype_code)
      data_blup$ID = rownames(data_blup)
      data_blup$y=data_blup$`(Intercept)`
      data_blup=data_blup[,c('ID','y')]
      i=intersect(data_blup$ID,rownames(X))
      X = X[i,]
      data_blup=data_blup[i,]
      subK=K[i,i]

      m4 = relmatLmer(y ~ 1 + unlist(X) + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
      se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      #rownames(se4)=founders
      names(se4)=c('value','se','tvalue')
      se4[2,]$value=se4$value[-1]+se4$value[1]
      se4$variable_f=factor(c(0,1),levels=c(0,1))

      png(sprintf('GridLMM/effect_sizes/founder_ES/%s_%s_F_and_H_BLUP_effect_sizes_lme4qtl.png',name,method),width=1000,height=800)
      print(ggplot(se4,aes(x=variable_f,y=value)) +
      geom_point() +
      geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
      theme(axis.text.x=element_text(size=10)) +
      xlab("Allele") + ylab("Effect Size") +
      labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
      dev.off()
    }
  }
  ses[[count]]=list(SE=se4,method=method,qtl_id=name,snp=snp,pos=pos)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/founder_ES/F_and_H_not_S_QTL_SEs.rds')
