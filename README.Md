# QTL detection in a MAGIC population using bi-allelic, parental, and ancestral haplotype variation

"Modeling allelic diversity of multiparent mapping populations affects detection of quantitative trait loci"

Full paper here: [Odell et al. (2022)](https://academic.oup.com/g3journal/article/12/3/jkac011/6509518)

Companion paper here: [Hudson et al. (2022)](https://academic.oup.com/g3journal/article/12/3/jkac013/6520465)

600K Affymetrix SNP Array genotype data can be found on [Figshare](https://figshare.com/articles/dataset/Biogemma_BALANCE_MAGIC_600K_SNP_Genotype_Data/14903262).

Phenotype and environment data can be found [here](https://figshare.com/s/5ee8337defdef63b04ce).

## Contents

- `scripts` : R scripts
- `slurm-script`: bash scripts, mostly for running R scripts in `scripts` on slurm
- 'genotypes': Founder and haplotype probabilities and physical and genetic maps used
- `GridLMM`: QTL, effect sizes, permutation significance threshold, and kinship matrices
- `ibd_segments`: RefinedIBD results and haplotype blocks
- `run_magicsim`: Files for generating simulated MAGIC population. See [magicsim](https://github.com/sarahodell/magicsim) for more information.


## Generating Haplotype probabilities

### 1. RefinedIBD
- Run `slurm-scripts/refined_ibd.sh`.
- Uses `genotypes/600K/Biogemma_Founders_600K_Genotypes_AGPv4_phased_biallelic_trimmed.vcf.gz` and `ibd_segments/plink_600K_genetic.map` as inputs.
- Writes out to `ibd_segments/refinedibd/600K/Biogemma_600K_Founders_RefinedIBD_chr{chr}.ibd.gz`
- The output will need to be unzipped before moving on to the next step.

### 2. Building haplotype blocks
- Run `scripts/build_ibd_segments.R`. It takes the argument, c, for the chromosome.
- Uses the files `genotypes/qtl2/startfiles/Biogemma_pmapc{c}.csv` and `ibd_segments/refinedibd/Biogemma_600K_Founders_RefinedIBD_chr{c}.ibd`.
- Writes out ibd blocks to `ibd_segments/refinedibd/600K/bg{c}_refined_ibd_blocks.txt` and the ibd graph to 'ibd_segments/refinedibd/bg{c}_ibd_graph.rds'.
- The minimum number of haplotypes for each chromosome is useful information for downstream analyses, and is writen out to `genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt`

### 3. Generating haplotoype probabilities
- Run `scripts/build_haplotypes.R`. It takes arguments c, for the chromosome, and h for the haplotype group (h through 16, where h is the minimum haplotype for that chromosome, as output above).

- Uses founder probabilities, a physical map, and the ibd blocks generated above. These are found `genotypes/probabilities/geno_probs/raw/bg{c}_genoprobs.rds`, `genotypes/qtl2/startfiles/Biogemma_pmapc{c}.csv`, and `ibd_segments/refinedibd/600K/bg{c}_refined_ibd_blocks.txt`, respectively.

- Writes out haplotypes probabilities to `genotypes/probabilities/haplotype_probs/RefinedIBD_600K/raw_bg{c}_refined_ibd_haplogroup{h}_probs.rds`

### 4. Filtering
- Run `scripts/filter_haplotypes.R`. It takes arguments c, for the chromosome, and h for the haplotype group (h through 16, where h is the minimum haplotype for that chromosome, as output above).
- Uses a physical map and the ibd segments and haplotype probabilities generated above. These are found `genotypes/qtl2/startfiles/Biogemma_pmapc{c}.csv`, `ibd_segments/refinedibd/600K/bg{c}_refined_ibd_blocks.txt`, and `genotypes/probabilities/haplotype_probs/RefinedIBD_600K/raw/bg{c}_refined_ibd_haplogroup{h}_probs.rds`, respectively.
- Writes out filtered haplotype probabilities based on an r^2 LD cutoff of 0.95 to `genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg{c}_filtered_haplogroup{h}_probs.rds`. Dropped markers and markers with low representation are also written out in the same directory.
