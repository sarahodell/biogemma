#!/env/usr/bin Rscript

library('data.table')
library('dplyr')
library('tibble')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])


founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


#wgs=fread(sprintf('biogemma/Biogemma_WGS_founder_alleles_chr%s.txt',c),data.table=F,header=F)
#names(wgs)=c('chr','pos','ref','alt1','alt2',founders)
#wgs = wgs%>% mutate(marker=paste0('S',chr,'_',pos))
#original=dim(wgs)[1]
#Drop rows where there are more than two alleles.
#wgs=wgs[wgs$alt2=='.',]
#biallelic=dim(wgs)[1]
#print(sprintf('Dropped %.0f sites that were multiallelic',original-biallelic))
#wgs=wgs[,c('marker','chr','pos','ref','alt1','alt2',founders)]

#for(s in founders){
#    d = fread(sprintf('biogemma/founders/%s.diff.sites',s),data.table=F)
#    d$founder=rep(s,dim(d)[1])
#    d=d[d$CHROM==c,]
#    d=d[d$N_DISCORD!=0,]
#    sites=match(d$POS,wgs$pos)
#    sites=subset(sites,is.na(sites)==F)
#    wgs[c(sites),c(s)]=NA
#}
#discord=dim(wgs)[1]
#print(sprintf("%.0f discordant sites dropped",biallelic-discord))
#len=dim(wgs)[1]
#binfinal=ifelse(wgs[,7:22]=='0/0',0,ifelse(wgs[,7:22]=='1/1',1,ifelse(wgs[7:22]=='2/2',2,NA)))
#binfinal=as.data.frame(binfinal)
#names(binfinal)=founders
#binfinal$marker=wgs$marker
#binfinal$pos=wgs$pos

#Drop rows with more than 12 founders with missing data
#binfinal = binfinal[rowSums(is.na(binfinal))<=12,]

#binfinal=binfinal[,c('marker','pos',founders)]
#print(sprintf('Keeping %.0f sites',dim(binfinal)[1]))

#fwrite(binfinal,file=sprintf('Biogemma_WGS_all_alleles_final_chr%s.txt',c),row.names=F,quote=F,sep='\t')
binfinal=fread(sprintf('Biogemma_WGS_all_alleles_final_chr%s.txt',c),data.table=F)

#Turn all NA values to zero for calculating allele probs.
freq=rowMeans(binfinal[,3:18],na.rm=T)
for(i in 3:18){
    ind=which(is.na(binfinal[,i]))
    binfinal[ind,i]=freq[ind]
}


options(scipen=999)

#for each chromosome
#read in files
pr=readRDS(sprintf('Biogemma_082318/600K_DHgenoprobs/bg%s_0823_genoprobs.rds',c))
#Read in physical map
pmap=fread(sprintf('Biogemma_082318/Biogemma_pmap_c%s.csv',c),data.table=F)
pmap$pos=pmap$pos*1e6

wgslen=dim(binfinal)[1]
sub=as.matrix(binfinal[,c(founders)])
pos=binfinal$pos
markers=binfinal$marker
remove(binfinal)
all_probs=c()
for( i in seq(1,344)){
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

all_probs=as.data.frame(all_probs)
names(all_probs)=markers
rownames(all_probs)=names(pr[[1]][,1,1])
all_probs=rownames_to_column(all_probs,"sample")
all_probs=all_probs[,c('sample',markers)
print("Writing allele probs to file")
fwrite(all_probs,sprintf('Biogemma_102218/Biogemma_WGS_chr%s_allele_probs.txt',c),row.names=F,sep='\t',quote=F)