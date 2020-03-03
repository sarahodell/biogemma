
# Filter Correlated Markers in Haplotype Probability

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])


cutoff=0.95

geno=readRDS(sprintf('bg%s_genoprobs_010319.rds',chr))

m_names=names(geno[[1]][,1,1])
dropped=list()
count=1
dropped[[count]]=list(marker=m_names[1],linked=c())

keep=c(1)
start=c()
for(k in 1:16){start=c(start,as.vector(geno[[1]][,k,1]))}
size=dim(geno[[1]])[3]
for(i in 1:size){
  ind_max = c()
  for(k in 1:16){ind_max=c(ind_max,as.vector(geno[[1]][,k,i]))}
  if(cor(start,ind_max)<cutoff){
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


saveRDS(dropped,sprintf('bg%s_dropped_markers_genoprobs.rds',chr))
saveRDS(filtered_geno,sprintf('bg%s_filtered_genotype_probs.rds',chr))
sprintf("Keeping %0.f markers on chromosome %s with correlation filter of %.2f",length(keep),chr,cutoff)