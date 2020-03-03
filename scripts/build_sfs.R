#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('dplyr')

all_freqs=c()

for(i in seq(1,10)){
      fgeno=fread(sprintf('../../Biogemma_foundergenos/Founder_genos_chr%.0f_121718.csv',i),data.table=F)
      freq=apply(fgeno,MARGIN=2,function(x) sum(x=="B"))
      freq=freq[2:length(freq)]
      for(i in seq(1,length(freq))){freq[i]=ifelse(freq[i]>8,(16-freq[i])/16,freq[i]/16)}
      all_freqs=c(all_freqs,freq)
}

sfs=as.data.frame(table(all_freqs),stringsAsFactors=F)
names(sfs)=c('maf','freq')
sfs$freq=as.integer(sfs$freq)

print(sfs)

png('600K_Founders_FoldedSFS.png',width=800,height=800)
print(ggplot(sfs,aes(x=maf,y=freq)) + geom_bar(stat="identity") + xlab("Minor Allele Frequency") + ylab("Count") + theme_classic() + ggtitle("MAGIC Founders Site Frequency Spectrum of 600K SNPs"))
dev.off()


all_freqs=c()

for(i in seq(1,10)){
      geno=fread(sprintf('../../Biogemma_DHgenos/DH_geno_chr%.0f_121718.csv',i),data.table=F)
      freq=apply(geno,MARGIN=2,function(x) sum(x=="B"))
      freq=freq[2:length(freq)]
      for(i in seq(1,length(freq))){freq[i]=ifelse(freq[i]>171.5,(343-freq[i])/343,freq[i]/343)}
      all_freqs=c(all_freqs,freq)
}


bins=seq(0,0.5,by=0.05)
freqbin=cut(all_freqs,bins)
sfs=as.data.frame(table(freqbin),stringsAsFactors=F)
names(sfs)=c('maf','freq')
sfs$freq=as.integer(sfs$freq)

print(sfs)

png('600K_DH_FoldedSFS.png',width=800,height=800)
print(ggplot(sfs,aes(x=maf,y=freq)) + geom_bar(stat="identity") + xlab("Minor Allele Frequency") + ylab("Count") + theme_classic() + ggtitle("MAGIC DH Site Frequency Spectrum of 600K SNPs"))
dev.off()