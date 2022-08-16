
# Filter Correlated Markers in Haplotype Probability

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])


library('data.table')

cutoff=0.95

genofile=fread(sprintf('DH_geno_chr%s_121718.csv',chr),data.table=F,stringsAsFactors=F)
print(dim(genofile))
inds=genofile$ind
variant=apply(genofile[,2:dim(genofile)[2]],MARGIN=2,function(x) length(unique(x))>1)
variant_names=names(variant[variant==T])
genofile=genofile[,colnames(genofile) %in% variant_names]
m_names=names(genofile)[2:dim(genofile)[2]]
genofile$ind=inds
genofile=genofile[,c('ind',m_names)]
print(dim(genofile))
fwrite(genofile,sprintf('DH_geno_chr%s_121718.csv',chr),row.names=F,quote=F,sep=',')
# Switch to 0 and 1 for reference and alternate alleles
geno=apply(genofile[,2:dim(genofile)[2]],MARGIN=2,function(x) ifelse(x=="A",0,1))
m_names=names(genofile)[2:dim(genofile)[2]]

#Filter out invariant sites
#mono=apply(geno,MARGIN=2,function(x) length(unique(x))>1)
#print(sum(mono))
m_names=m_names[mono]
geno=geno[,mono]
geno=as.data.frame(geno,stringsAsFactors=F)

geno$ind=genofile$ind
print(dim(geno))


n_ind=dim(geno)[1]
dropped=list()
count=1
dropped[[count]]=list(marker=m_names[1],linked=c())

keep=c(1)
start=c()
start=c(as.vector(geno[,2]))
size=dim(geno)[2]
for(i in 3:size){
  ind_max = c()
  ind_max=as.vector(geno[,i])
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


filtered_geno=geno[,keep]
print(dim(filtered_geno))

saveRDS(dropped,sprintf('bg%s_dropped_markers_600K.rds',chr))
fwrite(filtered_geno,sprintf('bg%s_filtered_600K.csv',chr),row.names=F,quote=F)
sprintf("Keeping %0.f markers on chromosome %s with correlation filter of %.2f",length(keep),chr,cutoff)
