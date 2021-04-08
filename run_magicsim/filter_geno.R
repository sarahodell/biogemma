
# Filter Correlated Markers in Haplotype Probability

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
rep=as.numeric(args[[2]])

library('data.table')

#r2 cutoff of 0.95
cutoff=0.95
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra",
           "OH43_inra","A654_inra","FV2_inra","C103_inra",
           "EP1_inra","D105_inra","W117_inra","B96","DK63",
           "F492","ND245","VA85")
geno=readRDS(sprintf('qtl2_files/MAGIC_DHSim_rep%.0f_c%s_genoprobs.rds',rep,chr))
dimnames(geno[[1]])[[2]]=founders
pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
pmap=pmap[order(pmap$pos),]
rownames(pmap)=seq(1,dim(pmap)[1])
pmap=pmap[pmap$marker %in% dimnames(geno[[1]])[[3]],]

# Drop lines that don't have both phenotype and genotype data
#K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
#rownames(K)=K[,1]
#rownames(K)=gsub("-",".",rownames(K))
#K=as.matrix(K[,-1])
#colnames(K)=rownames(K)

#phenotypes=fread('../GridLMM/phenotypes_asi.csv',data.table=F)
#phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
#pheno_IDs=unique(phenotypes$Genotype_code)

geno=lapply(geno,function(x) x[1:325,,])


#geno = lapply(geno,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.95,1,ifelse(x[,i]<=0.05,0,x[,i]))))
#for(i in 1:16){dimnames(geno[[i]])[[2]]=dimnames(geno[[i]])[[2]]}

geno=lapply(geno, function(x) x[,,pmap$marker])
#geno=sapply(geno,function(x) x[,,])
size=dim(geno[[1]])[3]
#Remove sites with low founder representation
# drop sites with summed founder prob of less than 1
f_sums=sapply(seq(1,size),function(x) colSums(geno[[1]][,,x]))
low_rep=which(f_sums<1,arr.ind=T)
low_rep=as.data.frame(low_rep,stringsAsFactors=F)
which_f=unique(low_rep$row)

dropped=vector("list",length=16)
for(f in which_f){
  cols=low_rep[low_rep$row==f,]$col
  #geno[[1]][,f,cols]=NA
  dropped[[f]]=dimnames(geno[[1]])[[3]][cols]
}

# drop sites with less than 5 lines with prob greater than .8
size=dim(geno[[1]])[3]
f_sums2=sapply(seq(1,size),function(x) colSums(geno[[1]][,,x]>=0.8))
low_rep2=which(f_sums2<5,arr.ind=T)
low_rep2=as.data.frame(low_rep2,stringsAsFactors=F)
which_f2=unique(low_rep2$row)
for(f in which_f2){
  cols=low_rep2[low_rep2$row==f,]$col
  #geno[[1]][,f,cols]=NA
  dropped[[f]]=unique(c(dropped[[f]],dimnames(geno[[1]])[[3]][cols]))
}

saveRDS(dropped,sprintf('qtl2_files/dropped/bg%s_low_rep_markers_rep%.0f.rds',chr,rep))

size=dim(geno[[1]])[3]
m_names=names(geno[[1]][1,1,])
dropped=list()
count=1
dropped[[count]]=list(marker=m_names[1],linked=c())

keep=c(1)
start=c()
for(k in 1:16){start=c(start,as.vector(geno[[1]][,k,1]))}

for(i in 1:size){
  ind_max = c()
  for(k in 1:16){ind_max=c(ind_max,as.vector(geno[[1]][,k,i]))}
  if((cor(start,ind_max,use="complete.obs")**2)<cutoff){
    keep=c(keep,i)
    start=ind_max
    start_ind=i
    count=count+1
    dropped[[count]]=list(marker=c(m_names[start_ind]),linked=c())
  }
  else{
     dropped[[count]]$linked=c(dropped[[count]]$linked,m_names[i])
  }
}


filtered_geno=list()
for(i in 1:16){filtered_geno[[i]]=geno[[1]][,i,keep]}
names(filtered_geno)=founders

saveRDS(dropped,sprintf('qtl2_files/dropped/bg%s_rep%.0f_dropped_markers.rds',chr,rep))
saveRDS(filtered_geno,sprintf('qtl2_files/filtered/bg%s_rep%.0f_filtered_genotype_probs.rds',chr,rep))
sprintf("Keeping %0.f markers on chromosome %s with correlation filter of %.2f",length(keep),chr,cutoff)
