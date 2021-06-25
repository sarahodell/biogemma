###Make Rqtl2 control file

library("qtl2")
library('data.table')
# Patho to sample genotype file
#geno_file_path="qtl2_files"

# Path to founder genotype file
fgeno_file_path="../../genotypes/qtl2/Biogemma_foundergenos"

# Cross type (see manual)
crosstype="riself16"


map_path="../../genotypes/qtl2/startfiles"
#Covariate file
covar=NULL

#Info on crossing scheme (see manual)
#crossinfo="Sim_cross_info.csv"

# Coding of genotypes (A for homozygous ref alleles, 2 for het, 3 for homozygous alt allele)
genocodes=c(A=1L,B=3L)
cross_info=fread('founder_cross_info.txt',data.table=F)
#Names of the founder lines
alleles=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

for(r in 1:100){
  alleles=unname(unlist(cross_info[cross_info$V1==r,2:17]))
  cross_file=data.frame(paste0('Sim',seq(1,344)),1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,stringsAsFactors=F)
  names(cross_file)=c('ind',alleles)
  crossinfo=sprintf("Sim%.0f_cross_info.csv",r)
  fwrite(cross_file,crossinfo,row.names=F,quote=F,sep=',')

  for (i in 1:10){
    #Name of control file
    control_file=sprintf("qtl2_files/MAGIC_DHSim_rep%.0f_c%.0f.json",r,i)

    # Name of sample genotype file
    geno_filename=sprintf("MAGIC_DHSimAll_rep%.0f_chr%.0f.csv",r,i)
    # Name of founder genotype file
    fgeno_filename=sprintf("Founder_genos_chr%.0f.csv",i)
    # Name of physical map
    pmap=sprintf("Biogemma_pmap_c%.0f.csv",i)
    # Name of genetic map
    gmap=sprintf("Biogemma_gmap_c%.0f.csv",i)


    # Brief description of file
    description=sprintf("400 Simulated DH MAGIC Lines from 16 Founders, Chromosome %.0f",i)

    write_control_file(output_file =control_file,
    crosstype = crosstype,geno_file=geno_filename,
    founder_geno_file = sprintf("%s/%s",fgeno_file_path,fgeno_filename),
    gmap_file=sprintf("%s/%s",map_path,gmap),pmap_file = sprintf("%s/%s",map_path,pmap),covar_file=covar,crossinfo_file = crossinfo,
    geno_codes=genocodes,
    alleles=alleles,
    sep=",",na.strings=c("NA"),comment.char="#",geno_transposed = FALSE,
    founder_geno_transposed = FALSE,description=description
    )
  }
}
