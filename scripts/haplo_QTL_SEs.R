#!/usr/bin/env Rscript


library('data.table')
library('lme4')
library('lme4qtl')
library('GridLMM')
library('ggplot2')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env=paste0(qtl$Phenotype,'_',qtl$Environment)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="Haplotype_probs",]
#qtl=qtl[qtl$Phenotype==pheno & qtl$Environment==env & qtl$Chromosome==as.numeric(chr),]
rownames(qtl)=seq(1,nrow(qtl))
ses=list()
count=1
for(i in 1:nrow(qtl)){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
  results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=results[results$SNP==snp,]$HAPGRP
  haplo_probs = readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
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
    m0 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
    se4=as.data.frame(summary(m0)$coef,stringsAsFactors=F)
    names(se4)=c('value','se','tvalue')
    rownames(se4)=paste0("HAPGRP_",seq(1,h))
    se4$hapgrp=rownames(se4)
    se4$variable_f=factor(se4$hapgrp,levels=paste0("HAPGRP_",seq(1,h)))

    null_model = GridLMM_ML(y~1 + (1|Genotype_code),phenotype,relmat=list(Genotype_code=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m0))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'Genotype_code'

    Y=as.matrix(phenotype$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    corr=cor(betas,se4$value)

  }
  else{
    phenotype = phenotype[!is.na(phenotype$y),]
    #phenotype$y=phenotype$y-mean(phenotype$y)
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    i=intersect(data_blup$ID,rownames(X))
    X = X[i,]
    data_blup=data_blup[i,]
    subK=K[i,i]
    m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK))

    se4=as.data.frame(summary(m0)$coef,stringsAsFactors=F)
    rownames(se4)=paste0("HAPGRP_",seq(1,h))
    names(se4)=c('value','se','tvalue')
    se4$hapgrp=rownames(se4)
    se4$variable_f=factor(se4$hapgrp,levels=paste0("HAPGRP_",seq(1,h)))
    #se4=se4[order(se4$value),]
    #GRIDLMM
    null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m0))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'ID'


    Y=as.matrix(data_blup$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    corr=cor(betas,se4$value)
  }
  ses[[count]]=list(SE=se4,qtl_id=name,snp=snp,pos=pos,corr=corr)
  count=count+1
    # pheno = cbind(pheno[match(rownames(X),pheno$Genotype_code),1:3],X_diff_mean)
    # model with no intercept
}

saveRDS(ses,'GridLMM/effect_sizes/Haplotype_prob_QTL_SEs.rds')

# model with intercept, dropping founder 1
#X_diff_mean = X %*% contr.sum(ncol(X))
#colnames(X_diff_mean) = paste0('F',1:ncol(X_diff_mean))
#m1 = relmatLmer(y ~ X_diff_mean + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
#coefs1 = data.frame(summary(m1)$coef)[-1,]

# estimating SE for F16
#coefs1$Precision = 1/coefs1[,2]^2
#coefs1$nobs = colSums(X)[-16]
#m1_slope = lm(Precision~nobs,coefs1[-1,])
#with(coefs1,plot(Precision~nobs))
#abline(m1_slope)
#F16 = -sum(coefs1[,1])
#F16_SE = 1/sqrt(predict(m1_slope,newdata = list(nobs = colSums(X)[16])))

# model with intercept, dropping founders 1 and 16
#X_drop1_diff_mean = X[,-1] %*% contr.sum(ncol(X[,-1]))
#colnames(X_drop1_diff_mean) = paste0('F',1+1:ncol(X_drop1_diff_mean))
#m2 = relmatLmer(male_flowering_d6 ~ X_drop1_diff_mean + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
#summary(m2)$coef
#coefs2 = data.frame(summary(m2)$coef)[-1,]

# estimating SE for F16
#coefs2$Precision = 1/coefs2[,2]^2
#coefs2$nobs = colSums(X)[-c(1,16)]
#m2_slope = lm(Precision~nobs,coefs2)
#with(coefs2,plot(Precision~nobs))
#abline(m2_slope)
#F16 = -sum(coefs2[,1])
#F16_SE = 1/sqrt(predict(m2_slope,newdata = list(nobs = colSums(X)[2])))

#fwrite(coefs,sprintf('GridLMM/effect_sizes/'))

#ibd=fread(sprintf('ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',chr),data.table=F)
#pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
#pos=pmap[pmap$marker==snp,]$pos
#ibd_seg=unlist(unname(ibd[ibd$start<=pos & ibd$end>pos,founders]))
#h=max(ibd_seg)
#haps=paste0("HAPGRP_",seq(1,h))
#drop2=which(colSums(X)<=5)
#if(length(drop2)!=0){
#  X_filt=X[,-drop2]
#}
#hmax=apply(X_filt,MARGIN=1,which.max)
#haps=haps[-drop2]
#data_blup=subset(data_blup,ID %in% rownames(X))

#data_blup$haplotype=haps[hmax]

#data_blup$haplotype_f=factor(data_blup$haplotype,levels=haps)
#m4 = relmatLmer(y ~ 0 + haplotype_f + (1|ID),data=data_blup,relmat = list(ID=K))
#se4=as.data.frame(summary(m4)$coef)
#names(se4)=c('value','se','tvalue')
#se4=se4[order(se4$value),]
#se4=as.data.frame(emmeans(m4,'haplotype_f'),stringsAsFactors=F)
#se3$founder=sort(f2)
#se4=se4[order(se4$emmean),]
#rownames(se4)=seq(1,nrow(se4))
#se4$variable=sapply(seq(1,nrow(se4)),function(x) strsplit(rownames(se4)[x],'founder')[[1]][2])
#se4$founder=as.character(se4$founder)
#se4$variable_f=factor(se4$founder,levels=se4$founder)

#founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#hapseq=seq(1,h)#[-drop2]
#founder_name=c(rep(founders,length(hapseq)))
#haplos=c()

#for(i in hapseq){haplos=c(haplos,rep(i,16))}
#founder_prop=c()
#for(j in hapseq){#
#  members=founders[which(ibd_seg==j)]
  #sub=se4[se4$haplotype==paste0('HAPGRP_',j),]
  #f=unique(se4$haplotype)
#  no=length(members)
#  for(k in founders){
#    if(k %in% members){
#      founder_prop=c(founder_prop,1/no)
#    }
#    else{
#      founder_prop=c(founder_prop,0)
#    }
#  }
#}
#fprop=data.frame(founder=founder_name,hapgrp=haplos,founder_prop=founder_prop,stringsAsFactors=F)

#h_value=c()
#h_low=c()
#h_high=c()
#for(o in hapseq){
#  v=se4[se4$hapgrp==paste0("HAPGRP_",o),]$value
#  low=v - (2*se4[se4$hapgrp==paste0("HAPGRP_",o),]$se)
#  high=v + (2*se4[se4$hapgrp==paste0("HAPGRP_",o),]$se)
#  h_value=c(h_value,rep(v,16))
#  h_low=c(h_low,rep(low,16))
#  h_high=c(h_high,rep(high,16))
#}
#fprop$h_value=h_value
#fprop$lower.CL=h_low
#fprop$higher.CL=h_high
#fprop$h_value_adj=fprop$founder_prop * fprop$h_value
#fprop$founder_f=factor(fprop$founder,levels=founders)

#colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
#rownames(colorcodes)=colorcodes$founder
#colorcodes=colorcodes[founders,]

#png('GridLMM/effect_sizes/qDTA3_1_ALL_haplotype_effect_sizes_lme4qtl.png',width=1000,height=800)
#print(ggplot(fprop,aes(x=as.factor(hapgrp),y=h_value_adj,fill=founder_f)) +
# geom_bar(position="stack",stat="identity")+
#  geom_errorbar(aes(ymin=lower.CL,ymax=higher.CL)) +
#  scale_fill_manual(values=colorcodes[levels(fprop$founder_f),]$hex_color,labels=levels(fprop$founder_f))+
#  ylab("Haplotype Effect Size") +
#  xlab("Haplotype Groups") + ggtitle("Haplotype BLUP Effect Sizes of qDTA3_1") +
#   theme(axis.text.x=element_text(size=8)))
#dev.off()


#png('GridLMM/effect_sizes/qDTA3_1_BLOIS_2017_OPT_effect_sizes_lme4qtl.png',width=1000,height=800)
#print(ggplot(se4,aes(x=variable_f,y=emmean)) +
#geom_point() +
 #geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) +
#theme(axis.text.x=element_text(size=10)) +
#xlab("Founder") + ylab("Effect Size (ggd)") +
#labs(title="Days to Anthesis, BLOIS_2017_OPT Effect Size Estimates",subtitle="qDTA3_1"))
#dev.off()
