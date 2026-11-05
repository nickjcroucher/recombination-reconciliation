# Pangenome structure and the homologous recombination / prophage conflict

Reproducible companion analysis testing whether the intragenomic conflict between
homologous recombination (HR) and prophage, described by
[Kwun & Croucher (2026)](https://doi.org/10.1093/ismejo/wrag256), can explain the
pangenome-scale variation across bacterial species reported by
[Horsfield et al. (2026)](https://doi.org/10.64898/2026.02.04.703695).

**[Render the full report](pangenome_recombination_conflict.Rmd)** - everything is
computed live from the three raw CSVs in `data/raw/`; no merged/derived data files
are committed to this repository.

## Structure

```
.
├── pangenome_recombination_conflict.Rmd   Main analysis document
├── R/
│   ├── data_prep.R                        Taxon matching / dataset construction
│   └── helpers.R                          Palettes, theme, plotting functions
├── data/raw/                              The three raw input CSVs (see data/raw/README.md)
└── .github/workflows/render.yml           CI: knits the Rmd on every push/PR
```

## Key findings

Flow of the analysis (figure numbers refer to the main document):

1. **Pangenome size tracks prophage density** (Figure 1; *r* = 0.85). Rare-gene
   count shows the same pattern, and transformable taxa have smaller pangenomes
   on average, but only because they also carry fewer prophage (Supplementary
   Figures S1-S2).
2. **The proportion of fast genes and the raw Y-intercept of the core/accessory
   distance curve are *not* explained by prophage density alone** (Supplementary
   Figures S3-S4) - but genome-size-scaling the Y-intercept (converting a per-site
   rate into an expected number of rapidly-diverging sites genome-wide) **does**
   recover a significant relationship with prophage density (Figure 2A; *r* = 0.64).
3. **The proportion of fast genes is instead explained by HR probability**
   (Figure 2B; *r* = -0.68), in the opposite direction to what would be expected if
   HR and prophage accumulation were the same underlying process.
4. **Basal gene turnover is unrelated to prophage density** (Supplementary Figure
   S5) **but rises significantly with HR probability** (Figure 2C; *r* = 0.66),
   consistent with transformation being an important mechanism facilitating the
   exchange of accessory loci that are not themselves autonomously mobile.

## Scope note: what counts as "prophage" here

Prophage density (`lambda_hat`) is restricted throughout to geNomad-detected,
hallmark-gene-rich *Caudoviricetes* provirus regions - large, plausibly-intact
temperate phage capable of lytic induction and host-cell death, matching the scope
of Kwun & Croucher's own cross-taxon comparison. Satellite phage/phage-related
chromosomal islands (PRCIs), *Inoviridae*, and plasmids are excluded; see
`data/raw/README.md` for details.

## Uncertainty / error bars

Every correlation plot shows 95% CI error bars on any axis for which the source
data report one: `lambda_hat` (Poisson fit, Kwun & Croucher), `Recombination
probability` (GLMM asymptotic CI, Kwun & Croucher), and `rate_genes1_median`,
`prop_genes2_median`, and the Supplementary File 8 Y-intercept (PopPUNK-mod
credible intervals / curve-fit CIs, Horsfield et al.). `pan_genes`, `rare_genes`,
and `genome_size` are reported by Horsfield et al. as point estimates only, so no
error bars are fabricated for these - this is noted explicitly in each affected
figure caption.

## Reproducing locally

```r
install.packages(c("tidyverse","ggrepel","cowplot","broom","knitr","kableExtra","scales","rmarkdown"))
rmarkdown::render("pangenome_recombination_conflict.Rmd")
```

## License

GNU GENERAL PUBLIC LICENSE Version 3 - see `LICENSE`.
