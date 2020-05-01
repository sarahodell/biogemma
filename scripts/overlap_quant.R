#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
c=as.character(args[[1]])

library('data.table')

overlap=fread(sprintf('ibd_segments/comparison/c%s_germline_overlap.bed',c),data.table=F)
names(overlap)=c('chr','start','end','name')
overlap$start=as.numeric(overlap$start)
overlap$end=as.numeric(overlap$end)
# Total IBD region of chromosome
pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)
chrom_ibd_perc=round((sum(overlap$end - overlap$start)/max(pmap$pos)) * 100,2)
total_chrom_ibd=sum(overlap$end - overlap$start)

# GERMLINE % Overlap with 600K IBD
sibd=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_600K_germline_IBD_chr%s.bed',c),data.table=F)
names(sibd)=c('chr','start','end','name')
sibd$start=as.numeric(sibd$start)
sibd$end=as.numeric(sibd$end)
overlap_600k=round((sum(overlap$end - overlap$start)/sum(sibd$end - sibd$start)) * 100,2)

#600K Summary
chrom_600k_perc=round((sum(sibd$end - sibd$start)/max(pmap$pos)) * 100,2)
chrom_600k_ibd=sum(sibd$end - sibd$start)

# GERMLINE % Overlap with WGS IBD
wibd=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_germline_IBD_chr%s.bed',c),data.table=F)
names(wibd)=c('chr','start','end','name')
wibd$start=as.numeric(wibd$start)
wibd$end=as.numeric(wibd$end)
overlap_wgs=round((sum(overlap$end - overlap$start)/sum(wibd$end - wibd$start)) * 100,2)

#WGS Summary
chrom_wgs_perc=round((sum(wibd$end - wibd$start)/max(pmap$pos)) * 100,2)
chrom_wgs_ibd=sum(wibd$end - wibd$start)


# IBDSeq Summary
ibdseq=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_IBDSeq_chr%s.bed',c),data.table=F)
names(ibdseq)=c('chr','start','end','name')
ibdseq$start=as.numeric(ibdseq$start)
ibdseq$end=as.numeric(ibdseq$end)
chrom_ibdseq_ibd=sum(ibdseq$end - ibdseq$start)

# RefinedIBD Summary
refined=fread(sprintf('ibd_segments/comparison/bedfiles/Biogemma_Founders_WGS_RefinedIBD_chr%s.bed',c),data.table=F)
names(refined)=c('chr','start','end','name')
refined$start=as.numeric(refined$start)
refined$end=as.numeric(refined$end)
chrom_refined_ibd=sum(refined$end - refined$start)

# GERMLINE 600K Overlap with IBDSeq
ibdseq_600K_overlap=fread(sprintf('ibd_segments/comparison/c%s_600K_ibdseq_overlap.bed',c),data.table=F)
names(ibdseq_600K_overlap)=c('chr','start','end','name')
ibdseq_600K_overlap$start=as.numeric(ibdseq_600K_overlap$start)
ibdseq_600K_overlap$end=as.numeric(ibdseq_600K_overlap$end)
ibdseq_600K_bp = sum(ibdseq_600K_overlap$end - ibdseq_600K_overlap$start)

overlap_1=round((sum(ibdseq_600K_overlap$end - ibdseq_600K_overlap$start)/sum(sibd$end - sibd$start)) * 100,2)
overlap_2=round((sum(ibdseq_600K_overlap$end - ibdseq_600K_overlap$start)/sum(ibdseq$end - ibdseq$start)) * 100,2)


# GERMLINE 600K Overlap with RefinedIBD
refined_600K_overlap=fread(sprintf('ibd_segments/comparison/c%s_germline_600K_RefinedIBD_overlap.bed',c),data.table=F)
names(refined_600K_overlap)=c('chr','start','end','name')
refined_600K_overlap$start=as.numeric(refined_600K_overlap$start)
refined_600K_overlap$end=as.numeric(refined_600K_overlap$end)
refined_600K_bp = sum(refined_600K_overlap$end - refined_600K_overlap$start)

overlap_5=round((sum(refined_600K_overlap$end - refined_600K_overlap$start)/sum(sibd$end - sibd$start)) * 100,2)
overlap_6=round((sum(refined_600K_overlap$end - refined_600K_overlap$start)/sum(refined$end - refined$start)) * 100,2)

# GERMLINE WGS Overlap with IBDSeq
ibdseq_wgs_overlap=fread(sprintf('ibd_segments/comparison/c%s_germline_wgs_ibdseq_overlap.bed',c),data.table=F)
names(ibdseq_wgs_overlap)=c('chr','start','end','name')
ibdseq_wgs_overlap$start=as.numeric(ibdseq_wgs_overlap$start)
ibdseq_wgs_overlap$end=as.numeric(ibdseq_wgs_overlap$end)
ibdseq_wgs_bp = sum(ibdseq_wgs_overlap$end - ibdseq_wgs_overlap$start)

overlap_3=round((sum(ibdseq_wgs_overlap$end - ibdseq_wgs_overlap$start)/sum(wibd$end - wibd$start)) * 100,2)
overlap_4=round((sum(ibdseq_wgs_overlap$end - ibdseq_wgs_overlap$start)/sum(ibdseq$end - ibdseq$start)) * 100,2)

header=sprintf('Chromosome %s IBD Comparison Summary\n\n',c)
summary=paste0(header,'Total GERMLINE 600K vs. GERMLINE WGS Consensus (bp): ',total_chrom_ibd,
'\nPercentage GERMLINE 600K Overlap: ',overlap_600k,
'\nPercentage Consensus Overlap with GERMLINE WGS: ',overlap_wgs,
'\n\n Total GERMLINE 600K vs. IBDSeq WGS Consensus (bp): ',ibdseq_600K_bp,
'\nPercentage Consensus Overlap with IBDSeq WGS: ',overlap_2,
'\nPercentage Consensus Overlap with GERMLINE 600K: ',overlap_1,
'\n\nTotal GERMLINE WGS vs. IBDSeq WGS Consensus (bp): ',ibdseq_wgs_bp,
'\nPercentage Consensus Overlap with IBDSeq WGS: ',overlap_4,
'\nPercentage Consensus Overlap with GERMLINE WGS: ',overlap_3,
'\n\nTotal GERMLINE 600K vs. RefinedIBD WGS Consensus (bp): ',refined_600K_bp,
'\nPercentage Consensus Overlap with RefinedIBD WGS: ',overlap_6,
'\nPercentage Consensus Overlap with GERMLINE 600K: ',overlap_5,
'\n\nTotal GERMLINE 600K IBD bp: ',chrom_600k_ibd,
'\n\nTotal GERMLINE WGS IBD bp: ',chrom_wgs_ibd,
'\n\nTotal IBDSeq WGS IBD bp: ', chrom_ibdseq_ibd,
'\n\nTotal RefinedIBD WGS bp: ', chrom_refined_ibd)
# 'Percentage of Chromosome in IBD: ',chrom_ibd_perc,'\n',
# '600K Percentage of Chromosome in IBD: ',chrom_600k_perc,
#,'\nWGS Percentage of Chromosome in IBD: ',chrom_wgs_perc,'\n')

fileConn<-file(sprintf("ibd_segments/comparison/c%s_summmary.txt",c))
writeLines(summary, fileConn)
close(fileConn)
