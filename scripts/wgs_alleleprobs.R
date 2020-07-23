#!/env/usr/bin Rscript
## HAHAHAHAHAH
library('data.table')
library('dplyr')
library('tibble')

args=commandArgs(trailingOnly=TRUE)
chr=as.character(args[[1]])

genofile=sprintf('genotypes/probabilities/geno_probs/raw/bg%s_genoprobs_010319.rds',chr)
#imputefile=sprintf('gzip -dc genotypes/WGS/Biogemma_WGS_all_alleles_final_chr%s.txt.gz',chr)
imputefile=sprintf('gzip -dc /scratch/sodell/Biogemma_WGS_all_alleles_final_chr%s.txt.gz',chr)
pmapfile=sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr)
outfile=sprintf('genotypes/probabilities/allele_probs/bg%s_wgs_alleleprobs_bimbam.txt',chr)

print("Starting allele probability imputation step")

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

binfinal=fread(cmd=imputefile,data.table=F,header=F)
names(binfinal)=c("marker","pos","alt1","ref",founders)
binfinal[,5:20] = ifelse(binfinal[,5:20]=="0/0",0,ifelse(binfinal[,5:20]=="1/1",1,NA))
binfinal$marker=paste0('S',chr,'_',binfinal$pos)
#Turn all NA values to zero for calculating allele frequency.
freq=rowMeans(binfinal[,5:20],na.rm=T)
print("Removing monomorphic SNPs")
binfinal=binfinal[freq>0 & freq<1,]
binfinal=as.data.frame(binfinal,stringsAsFactors=F)
rownames(binfinal)=seq(1,dim(binfinal)[1])

for(i in 5:20){
    ind=which(is.na(binfinal[,i]))
    binfinal[ind,i]=freq[ind]
}

options(scipen=999)

#for each chromosome
#read in files
pr=readRDS(genofile)
#Read in physical map
pmap=fread(pmapfile,data.table=F)

print("Ordering the samples based on the kinship matrix")
K=fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


print("Calculating allele probabilities")
wgslen=dim(binfinal)[1]
sub=as.matrix(binfinal[,c(founders)])
pos=binfinal$pos
markers=binfinal$marker
samples=rownames(K)
size=length(samples)
pr[[1]]=pr[[1]][samples,,]
#samples=names(pr[[1]][,1,1])
alt1=binfinal$alt1
ref=binfinal$ref
remove(binfinal)
all_probs=c()
for( i in seq(1,size)){
    print(i)
    ind=pr[[1]][i,,]
    ind=t(ind)
    ind=as.data.frame(ind)
    ind=rownames_to_column(ind,"marker")
    names(ind)=c("marker",founders)
    ind=merge(ind,pmap,by.x="marker",by.y="marker")
    ind=ind[order(ind$pos),]
    rownames(ind)=seq(1,nrow(ind))
    ind=ind[,c("marker","chr","pos",founders)]
    prlen=dim(ind)[1]
    f_probs=c()
    for(f in seq(1,16)){
        founder=founders[f]
        f_ind=ind[,c('pos',founder)]
        f_interp=approxfun(f_ind$pos,f_ind[,c(founder)],method="linear",yleft=unlist(unname(f_ind[1,c(founder)])),yright=unlist(unname(f_ind[prlen,c(founder)])))
        f_probs=rbind(f_probs,f_interp(pos))
    }
    f_probs=t(as.matrix(f_probs))
    allele_probs=rowSums(sub*f_probs)
    all_probs=rbind(all_probs,allele_probs)
}
print("Finished calculating allele probabilities")
all_probs_t=t(all_probs)
#remove(all_probs)
all_probs_t=as.data.frame(all_probs_t,stringsAsFactors=F)
names(all_probs_t)=samples
all_probs_t$marker=markers
#all_probs_t=rownames_to_column(all_probs_t,"marker")
all_probs_t$alt1=alt1
all_probs_t$ref=ref
all_probs_t=all_probs_t[,c('marker','alt1','ref',samples)]


print("Converting from alternate allele to minor allele")
freq=rowSums(all_probs_t[4:dim(all_probs_t)[2]])/(dim(all_probs_t)[2]-3)
min_allele=freq<=0.5

is_major=which(min_allele==F)

bimbam=function(i,df){
  line=df[i,]
  tmp=df[i,'alt1']
  df[i,'alt1']=df[i,'ref']
  df[i,'ref']=tmp
  df[i,4:dim(df)[2]]=1-df[i,4:dim(df)[2]]
  return(df)
}

df_edit <- bimbam(is_major,all_probs_t)

# Set frequencies from 0..2
df_edit[,4:dim(df_edit)[2]] = df_edit[,4:dim(df_edit)[2]]*2
df_edit_o=na.omit(df_edit)

print("Writing allele probs to file")
fwrite(df_edit_o,outfile,row.names=F,sep=',',quote=F,col.names=F)
#system(sprintf("gzip %s",outfile))
