#!/env/usr/bin Rscript
#install.packages("data.table")

library('data.table')
library('dplyr')
library('tibble')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])


founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


wgs=fread(sprintf('biogemma/Biogemma_WGS_founder_alleles_chr%s.txt',c),data.table=F,header=F)
names(wgs)=c('chr','pos','ref','alt1','alt2',founders)
wgs = wgs%>% mutate(marker=paste0('S',chr,'_',pos))
wgs=wgs[wgs$alt2=='.',]
wgs=wgs[,c('marker','chr','pos','ref','alt1','alt2',founders)]

#array=fread(sprintf('biogemma/Biogemma_600K_founder_alleles_chr%s.txt',c),data.table=F,header=F)
#names(array)=c('chr','pos','ref','alt1','alt2',founders)
#array = array%>% mutate(marker=paste0('S',chr,'_',pos))

#w=match(array$marker,wgs$marker)
#final=subset(wgs,!(rownames(wgs) %in% w))
#final=rbind(final,array)
#final=final[order(final$pos),]
#rownames(final)=seq(1,nrow(final))
#final=final[,c('marker','chr','pos','ref','alt1','alt2',founders)]

#remove(wgs)
#remove(array)
print("Got here")
for(s in founders){
    d = fread(sprintf('biogemma/founders/%s.diff.sites',s),data.table=F)
    d$founder=rep(s,dim(d)[1])
    d=d[d$CHROM==c,]
    d=d[d$N_DISCORD!=0,]
    print(dim(d))
    sites=match(d$POS,wgs$pos)
    sites=subset(sites,is.na(sites)==F)
    wgs[c(sites),c(s)]=NA
}
print("Got here. Discordant sites dropped")
len=dim(wgs)[1]
binfinal=ifelse(wgs[,7:22]=='0/0',0,ifelse(wgs[,7:22]=='1/1',1,ifelse(wgs[7:22]=='2/2',2,NA)))
binfinal=as.data.frame(binfinal)
names(binfinal)=founders
binfinal$marker=wgs$marker

#Drop rows with more than 12 founders with missing data
binfinal = binfinal[rowSums(is.na(binfinal))<=12,]

binfinal = binfinal %>% mutate(pos=strsplit(marker,'_')[[1]][2])
binfinal$pos=as.numeric(binfinal$pos)
binfinal=binfinal[,c('marker','pos',founders)]
# Drop rows where there are more than two alleles
print(dim(binfinal))

fwrite(binfinal,file=sprintf('Biogemma_WGS_all_alleles_final_chr%s.txt',c),row.names=F,quote=F,sep='\t')

