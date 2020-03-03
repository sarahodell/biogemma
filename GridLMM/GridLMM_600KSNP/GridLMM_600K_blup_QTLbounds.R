#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('ggplot2')

# Read in Kinship Matrix
K=fread(sprintf('../K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('../phenotypes_asi.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

#map=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

geno=fread(sprintf('DH_geno_chr%s_121718.csv',chr),data.table=F,stringsAsFactors=F)
X=sapply(seq(1,dim(geno)[1]),function(x) ifelse(geno[x,2:dim(geno)[2]]=='A',0,1))
X=t(X)
print(dim(X))
X=as.data.frame(X)
rownames(X)=geno$ind
colnames(X)=colnames(geno)[2:dim(geno)[2]]

X=X[rownames(X) %in% data_blup$ID,]
X=as.matrix(X)
data=data[data$ID %in% rownames(X),]

gwas_null = GridLMM_GWAS(
                        formula = y~1 + (1|ID),
                        test_formula = ~1,
                        reduced_formula = ~0,
                        data = data_blup,
                        weights = NULL,
                        X = X,
                        X_ID = 'ID',
                        h2_start = NULL,
                        h2_step = 0.01,
                        max_steps = 100,
                        relmat = list(ID=K),
                        centerX = TRUE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'ML',
                        mc.cores = cores,
                        verbose = FALSE
)

#saveRDS(gwas,sprintf('models/chr%s_%s_x_ALL_600KSNP_ML.rds',chr,pheno))

threshtable=fread('../threshold_table.txt',data.table=F)
cutoff=threshtable[threshtable$phenotype==pheno & threshtable$method=="600K_SNP" & threshtable$environment=="ALL",]$threshold

candidate_marker=gwas_null$results[which.min(gwas_null$results$p_value_ML),]$X_ID
min_pval=gwas_null$results[which.min(gwas_null$results$p_value_ML),]$p_value_ML

names=colnames(X)
isit=names==candidate_marker

location=which(isit==T,isit)
data_blup$candidate = X[,candidate_marker]
X_minus=X[,-location]


gwas_marker = GridLMM_GWAS(
                        formula = y~candidate + (1|ID),
                        test_formula = ~1,
                        reduced_formula = ~0,
                        data = data_blup,
                        weights = NULL,
                        X = X_minus,
                        X_ID = 'ID',
                        h2_start = NULL,
                        h2_step = 0.01,
                        max_steps = 100,
                        relmat = list(ID=K),
                        centerX = TRUE,
                        scaleX = FALSE,
                        fillNAX = FALSE,
                        method = 'ML',
                        mc.cores = cores,
                        verbose = FALSE
)

gwas_null_drop=gwas_null$results[gwas_null$results$X_ID!=candidate_marker,]

pmap=fread(sprintf('../../qtl2_startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
p=match(gwas_null$results$X_ID,pmap$marker)
phy_pos=pmap[p,]$pos
chrom=data.frame(chrm=chr,ID=gwas_null$results$X_ID,pos=phy_pos,pvalues=gwas_null$results$p_value_ML,stringsAsFactors = F)

min_pos=chrom[chrom$ID==candidate_marker,]$pos

chrom$log10p=-log10(chrom$pvalues)
chrom$sig = chrom$log10p >= cutoff
rownames(chrom)=seq(1,dim(chrom)[1])

chrom=chrom[chrom$ID!=candidate_marker,]

chrom$lrt=(-log10(pchisq(gwas_marker$results$ML_logLik-gwas_null_drop$ML_logLik,1,lower.tail=F)))

chrom=chrom[!is.na(chrom$pvalues),]
chrom=chrom[!is.na(chrom$pos),]

fwrite(chrom,sprintf('bound_images/%s_ALL_600K_SNP_chr%s_LRT.txt',pheno,chr),row.names=F,quote=F)

png(sprintf('bound_images/%s_ALL_600K_SNP_chr%s_LRT.png',pheno,chr),width=800,height=800)
print(ggplot(chrom,aes(x=pos,y=lrt)) + geom_point(aes(color=sig),alpha=0.5) + scale_color_manual(breaks=chrom$sig,values=c("FALSE"="black","TRUE"="red")) + geom_hline(yintercept=cutoff,color="green") + geom_hline(yintercept=(-log10(min_pval)),color="blue") + guides(color=F) + xlab("Position (Mb)") + ylab("-log10(p-value)") + ggtitle(sprintf("%s BLUP 600K LRT Chr %s",pheno,chr)) + theme_classic())
dev.off()