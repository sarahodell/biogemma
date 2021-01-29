library('data.table')
library('lme4')
library('lme4qtl')
library('ggplot2')
library('reshape2')
library('emmeans')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="Founder_probs",]

poi1="AX-91768118"
pos1=126077002

poi2="AX-91772415"
pos2=150347882

chr="8"
pheno="male_flowering_d6"
env="ALL"
name="male_flowering_d6"

founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
phenotype = subset(phenotype,Genotype_code %in% rownames(K))

#poi1
X = do.call(cbind,lapply(founder_probs,function(x) x[,poi1]))
colnames(X) = founders
rownames(X) = dimnames(founder_probs[[1]])[[1]]
phenotype = phenotype[!is.na(phenotype$y),]
phenotype$y=phenotype$y-mean(phenotype$y)
m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
data_blup = as.data.frame(ranef(m0)$Genotype_code)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]
X = X[rownames(X) %in% data_blup$ID,]

X_clean=t(sapply(seq(1,nrow(X)),function(x) ifelse(X[x,]>0.85,1,0)))
rownames(X_clean)=rownames(X)
drop2=which(colSums(X_clean)<5)
drop=which(rowSums(X_clean)<1)
X_filt=X_clean[-drop,-drop2]

if(length(drop2)!=0){
  f2=founders[founders!=founders[drop2]]
}else{
  f2=founders
}
#Using discrete founder variables
fmax=apply(X_filt,MARGIN=1,which.max)
data_blup1=data_blup[rownames(X_filt),]
data_blup1$founder=f2[fmax]
# These MITE probabilities are wrong ughhhhh
miteprob=fread('GridLMM/mite_probabilities.txt',data.table=F)
names(miteprob)=c('ID','mite')
miteprob$mite=ifelse(miteprob$mite>=0.85,1,0)
rownames(miteprob)=miteprob$ID
data_blup1$mite=miteprob[rownames(data_blup1),]$mite

K1=K[rownames(K) %in% rownames(X_filt),colnames(K) %in% rownames(X_filt)]

m1 = relmatLmer(y ~ 0 + founder + (1|ID),data=data_blup1,relmat = list(ID=K1))
se1=as.data.frame(emmeans(m1,'founder'),stringsAsFactors=F)
se1=se1[order(se1$emmean),]
rownames(se1)=seq(1,nrow(se1))

se1$variable_f=factor(se1$founder,levels=se1$founder)
se1$mite=has_mite[match(se1$founder,founders)]

m2 = relmatLmer(y ~ 0 + mite + founder + (1|ID),data=data_blup1,relmat = list(ID=K1))
se2=as.data.frame(emmeans(m2,'founder'),stringsAsFactors=F)
se2=se2[order(se2$emmean),]
rownames(se2)=seq(1,nrow(se2))

se2$variable_f=factor(se2$founder,levels=se2$founder)
se2$mite=has_mite[match(se2$founder,founders)]

png(sprintf('GridLMM/effect_sizes/%s_BLUP_effect_sizes_lme4qtl.png',poi1),width=1000,height=800)
print(ggplot(se2,aes(x=variable_f,y=emmean,color=mite)) +
geom_point() +
 geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) +
theme(axis.text.x=element_text(size=10)) +
xlab("Founder") + ylab("Effect Size (ggd)") +
labs(title="Days to Anthesis BLUP Effect Size Estimates",subtitle=subtitle=sprintf("Chr8: %.0f",pos1),color="MITE Present"))
dev.off()

#poi2
X = do.call(cbind,lapply(founder_probs,function(x) x[,poi2]))
colnames(X) = founders
rownames(X) = dimnames(founder_probs[[1]])[[1]]

X_clean=t(sapply(seq(1,nrow(X)),function(x) ifelse(X[x,]>0.85,1,0)))
rownames(X_clean)=rownames(X)
drop2=which(colSums(X_clean)<5)
drop=which(rowSums(X_clean)<1)

if(length(drop2)!=0){
  X_filt=X_clean[-drop,-drop2]
  f2=founders[founders!=founders[drop2]]
}else{
  X_filt=X_clean[-drop,]
  f2=founders
}
#Using discrete founder variables
fmax=apply(X_filt,MARGIN=1,which.max)
data_blup1=data_blup[rownames(X_filt),]
data_blup1$founder=f2[fmax]
# These MITE probabilities are wrong ughhhhh
miteprob=fread('GridLMM/mite_probabilities.txt',data.table=F)
names(miteprob)=c('ID','mite')
miteprob$mite=ifelse(miteprob$mite>=0.85,1,0)
rownames(miteprob)=miteprob$ID
data_blup1$mite=miteprob[rownames(data_blup1),]$mite

K1=K[rownames(K) %in% rownames(X_filt),colnames(K) %in% rownames(X_filt)]

m1 = relmatLmer(y ~ 0 + founder + (1|ID),data=data_blup1,relmat = list(ID=K1))
se1=as.data.frame(emmeans(m1,'founder'),stringsAsFactors=F)
se1=se1[order(se1$emmean),]
rownames(se1)=seq(1,nrow(se1))

se1$variable_f=factor(se1$founder,levels=se1$founder)
se1$mite=has_mite[match(se1$founder,founders)]

m2 = relmatLmer(y ~ 0 + mite + founder + (1|ID),data=data_blup1,relmat = list(ID=K1))
se2=as.data.frame(emmeans(m2,'founder'),stringsAsFactors=F)
se2=se2[order(se2$emmean),]
rownames(se2)=seq(1,nrow(se2))

se2$variable_f=factor(se2$founder,levels=se2$founder)
se2$mite=has_mite[match(se2$founder,founders)]

png(sprintf('GridLMM/effect_sizes/%s_BLUP_effect_sizes_lme4qtl.png',poi2),width=1000,height=800)
print(ggplot(se2,aes(x=variable_f,y=emmean,color=mite)) +
geom_point() +
 geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) +
theme(axis.text.x=element_text(size=10)) +
xlab("Founder") + ylab("Effect Size (ggd)") +
labs(title="Days to Anthesis BLUP Effect Size Estimates",subtitle=sprintf("Chr8: %.0f",pos2),color="MITE Present"))
dev.off()
