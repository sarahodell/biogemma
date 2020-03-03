#!/env/usr/bin Rscript

library('data.table')
library('dplyr')
library('tibble')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])


founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


wgs=fread(sprintf('../biogemma/Biogemma_WGS_founder_alleles_chr%s.txt',c),data.table=F,header=F)
names(wgs)=c('chr','pos','ref','alt1','alt2',founders)
wgs = wgs%>% mutate(marker=paste0('S',chr,'_',pos))
original=dim(wgs)[1]
#Drop rows where there are more than two alleles.
wgs=wgs[wgs$alt2=='.',]
biallelic=dim(wgs)[1]
print(sprintf('Dropped %.0f sites that were multiallelic',original-biallelic))
wgs=wgs[,c('marker','chr','pos','ref','alt1','alt2',founders)]

for(s in founders){
    d = fread(sprintf('../biogemma/founders/%s.diff.sites',s),data.table=F)
    d$founder=rep(s,dim(d)[1])
    d=d[d$CHROM==c,]
    d=d[d$N_DISCORD!=0,]
    sites=match(d$POS,wgs$pos)
    sites=subset(sites,is.na(sites)==F)
    wgs[c(sites),c(s)]=NA
}
discord=dim(wgs)[1]
print(sprintf("%.0f discordant sites dropped",biallelic-discord))
len=dim(wgs)[1]
binfinal=ifelse(wgs[,7:22]=='0/0',0,ifelse(wgs[,7:22]=='1/1',1,ifelse(wgs[7:22]=='2/2',2,NA)))
binfinal=as.data.frame(binfinal)
names(binfinal)=founders
binfinal$marker=wgs$marker
binfinal$pos=wgs$pos
binfinal$ref=wgs$ref
binfinal$alt1=wgs$alt1

#Drop rows with more than 12 founders with missing data
binfinal = binfinal[rowSums(is.na(binfinal))<=12,]

binfinal=binfinal[,c('marker','pos','alt1','ref',founders)]
print(sprintf('Keeping %.0f sites',dim(binfinal)[1]))

fwrite(binfinal,file=sprintf('Biogemma_WGS_all_alleles_final_chr%s.txt',c),row.names=F,quote=F,sep='\t')