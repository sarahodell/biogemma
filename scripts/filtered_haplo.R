
# Filter Correlated Markers in Haplotype Probability

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
hap_n=args[[2]]
#date=format(Sys.time(), "%m%d%y")

hap=readRDS(sprintf('bg%s_haplotype_probs_010319.rds',c))
pmap=fread(sprintf('../../qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
pmap=pmap[order(pmap$pos),]
rownames(pmap)=seq(1,dim(pmap)[1])
for(h in 1:length(hap)){
  hap_markers=dimnames(hap)[[2]]
}
hap=lapply(hap,function(x) x[,,pmap$marker]

m_names=names(hap[[1]][1,])
dropped=list()
count=1
dropped[[count]]=list(marker=c(m_names[1],linked=c())

keep=c(1)
start=c()
for(h in 1:hap_n){start=c(start,as.vector(hap[[h]][,1]))}
for(i in 2:dim(hap[[1]])[2]){
  ind_max = c()
  ind_max = for(h in 1:hap_n){ind_max=c(ind_max,as.vector(hap[[h]][,i]))}
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

filtered_hap=c()
for(h in 1:hap_n){filtered_hap[[h]]=hap[[h]][,keep]}


saveRDS(dropped,sprintf('bg%s_dropped_markers.rds'c))
saveRDS(filtered_hap_all,sprintf('bg%s_filtered_haplotype_probs.rds',c,date))
sprintf("Keeping %0.f markers on chromosome %s with correlation filter of %.2f",length(keep),chr,cutoff)
