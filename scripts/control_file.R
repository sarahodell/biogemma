###Make Rqtl2 control file

library("qtl2")

# Patho to sample genotype file
geno_file_path="../Biogemma_DHgeno"

# Path to founder genotype file
fgeno_file_path="../Biogemma_foundergenos"

# Cross type (see manual)
crosstype="riself16"

#Covariate file
covar=NULL

#Info on crossing scheme (see manual)
crossinfo="Biogemma_cross_info.csv"

# Coding of genotypes (A for homozygous ref alleles, 2 for het, 3 for homozygous alt allele)
genocodes=c(A=1L,B=3L)

#Names of the founder lines
alleles=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

for (i in 1:10){
    #Name of control file
    control_file=sprintf("Biogemma_c%.0f.json",i)

    # Name of sample genotype file
    geno_filename=sprintf("DH_geno_chr%.0f_121718.csv",i)
    # Name of founder genotype file
    fgeno_filename=sprintf("Founder_genos_chr%.0f_121718.csv",i)
    # Name of physical map
    pmap=sprintf("Biogemma_pmap_c%.0f.csv",i)
    # Name of genetic map
    gmap=sprintf("Biogemma_gmap_c%.0f.csv",i)


    # Brief description of file
    description=sprintf("344 DH MAGIC Lines from 16 Biogemma Founders, Chromosome %.0f",i)

    write_control_file(output_file =control_file,
    crosstype = crosstype,geno_file=sprintf("%s/%s",geno_file_path,geno_filename),
    founder_geno_file = sprintf("%s/%s",fgeno_file_path,fgeno_filename),
    gmap_file=gmap,pmap_file = pmap,covar_file=covar,crossinfo_file = crossinfo,
    geno_codes=genocodes,
    alleles=alleles,
    sep=",",na.strings=c("NA"),comment.char="#",geno_transposed = FALSE,
    founder_geno_transposed = FALSE,description=description
    )
}