
# Filter Correlated Markers in Haplotype Probability

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])


geno=readRDS(sprintf('bg%s_genoprobs_010319.rds',chr))

regrp_geno=list()
for(i in 1:16){regrp_geno[[i]]=geno[[1]][,i,]}

saveRDS(regrp_geno,sprintf('bg%s_reformat_genotype_probs.rds',chr))
