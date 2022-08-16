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
