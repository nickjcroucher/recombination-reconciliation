# Pangenome structure, homologous recombination and prophage conflict

Reproducible analysis testing the hypothesis that the **intragenomic
conflict between homologous recombination (HR) and prophage**, described by
Kwun & Croucher (2026, *ISME Journal*,
[doi:10.1093/ismejo/wrag256](https://doi.org/10.1093/ismejo/wrag256)),
explains the taxon-scale pangenome structure variation reported by
Horsfield *et al.* (2026, bioRxiv,
[doi:10.64898/2026.02.04.703695](https://doi.org/10.64898/2026.02.04.703695))
from fitting the *Pansim* / *PopPUNK-mod* framework to over 400 bacterial
species.

This repository follows the plotting conventions, palettes and general
structure of
[nickjcroucher/recombination-conflict](https://github.com/nickjcroucher/recombination-conflict),
re-using its `niche_palette` and `transformability_palette` colour schemes,
plotmath-formatted taxon labels, `geom_label_repel`-annotated scatterplots,
and violin/jitter comparison plots with Wilcoxon rank-sum annotation.

**No pre-merged or derived dataset is committed to this repository.** The
three files in `data/raw/` are the original, unmodified source files; every
merge, taxon match, and aggregation step is performed live by
`R/data_prep.R` when the main R Markdown document is knitted, so the full
provenance of every number in the analysis is traceable and reproducible.

## Repository structure

```
.
├── pangenome_recombination_conflict.Rmd   # Main analysis document (knit this)
├── R/
│   ├── data_prep.R                        # Builds the analysis datasets from data/raw/
│   └── helpers.R                          # Palettes, theme and plotting functions
├── data/
│   ├── raw/
│   │   ├── supplementary_file1_408species.csv            (Horsfield et al., Supp. File 1)
│   │   ├── supplementary_file8_negexp_fits.csv            (Horsfield et al., Supp. File 8)
│   │   └── recombination_conflict_comparison_plot.csv     (Kwun & Croucher)
│   └── README.md                          # Full data dictionary and merge logic
├── figures/                               # Populated when the Rmd is knitted
└── .github/workflows/render.yml           # CI: knits the Rmd on every push/PR
```

## Source data

| File | Origin | Content |
|---|---|---|
| `data/raw/recombination_conflict_comparison_plot.csv` | Kwun & Croucher | Per-taxon (n = 23) geNomad Caudoviricetes prophage density (`lambda_hat` + 95% CI) and probability of diversification through HR (`Recombination probability` + 95% CI), `Transformability`, `Niche` |
| `data/raw/supplementary_file1_408species.csv` | Horsfield *et al.*, Supplementary File 1 | PopPUNK-mod pangenome-dynamics parameters and WTBCluster pangenome composition statistics for 408 species |
| `data/raw/supplementary_file8_negexp_fits.csv` | Horsfield *et al.*, Supplementary File 8 | Negative-exponential core-vs-accessory distance curve fits (Plateau, Rate ratio, Y-intercept) for 22 species/datasets |

See `data/README.md` for the full column dictionary and taxon-matching
logic (direct matches, genus-level means for "spp." taxa, and the
*Enterococcus faecium* / `Enterococcus_B faecium` renaming).

## Reproducing the analysis

```r
install.packages(c("rmarkdown", "knitr", "kableExtra", "tidyverse",
                   "ggrepel", "cowplot", "scales", "broom"))
rmarkdown::render("pangenome_recombination_conflict.Rmd")
```

Or, from the command line:

```sh
Rscript -e 'rmarkdown::render("pangenome_recombination_conflict.Rmd")'
```

Continuous integration (`.github/workflows/render.yml`) knits the document
on every push and pull request against `main`, so any change that breaks
the pipeline — including changes to the raw data or the merge logic in
`R/data_prep.R` — will be caught automatically. The rendered HTML report is
uploaded as a build artifact.

## Summary of findings

| Relationship | r | p |
|---|---|---|
| Pangenome size (`pan_genes`) ~ prophage density (`lambda_hat`) | 0.85 | 3.7×10⁻⁷ |
| Pangenome size ~ HR probability | −0.48 | 0.021 |
| Proportion of fast-exchanged genes ~ HR probability | −0.68 | 3.4×10⁻⁴ |
| Basal gene turnover rate ~ HR probability | 0.66 | 6.2×10⁻⁴ |
| Core mutation rate ~ HR probability | −0.40 | 0.059 |
| (Y-intercept × genome size) ~ prophage density | 0.64 | 0.014 |

Pangenome size and the number of rare/accessory genes track **prophage
density**, not HR. The partitioning of the accessory genome into
"fast"/"slow" turnover compartments, and the basal turnover rate itself,
instead track **HR probability**. Naturally transformable taxa have
significantly smaller pangenomes (Wilcoxon test), but a joint linear model
(`pan_genes ~ lambda_hat + Transformability`) shows that transformability
adds no explanatory power once prophage density is accounted for
(p ≈ 0.61) — i.e. the effect of transformability on pangenome size is
mediated entirely through its association with lower prophage burden,
consistent with chromosomal curing by transformation.

## Uncertainty / error bars

Every scatterplot of two continuous variables shows 95% confidence/credible
interval error bars (vertical and/or horizontal via `geom_errorbar()` /
`geom_errorbarh()`) wherever the source data report one, via an explicit
lookup table in `R/helpers.R` (`ci_lookup`):

| Variable | CI source |
|---|---|
| `lambda_hat` | Poisson MLE 95% CI (`lambda_lwr`/`lambda_upr`) |
| `Recombination probability` | Asymptotic 95% CI from the GLMM (`asymp.LCL`/`asymp.UCL`) |
| `rate_genes1_median`, `core_mu_median`, `prop_genes2_median` | PopPUNK-mod 95% Bayesian credible intervals |

Variables with **no reported uncertainty** in the underlying supplementary
data — raw pangenome gene counts (`pan_genes`, `rare_genes`), `genome_size`,
and the negative-exponential curve-fit point estimates (`Plateau (i)`,
`Rate ratio (k)`, `Y-intercept (j)`) — are plotted without fabricated error
bars; each affected figure in the Rmd carries an explicit caption note
stating this rather than silently omitting the bars. Violin/jitter group
comparisons (`plot_violin_comparison()`) show the full distribution of
taxon-level values plus an annotated Wilcoxon test, which is the
uncertainty representation used for categorical comparisons in the source
`recombination-conflict` repository.

## Plotting conventions

- **One relationship per plot.** `Transformability` (colour) and `Niche`
  (shape) are used purely to annotate individual points; they are never
  used to sub-group the fitted regression. In `plot_correlation()`, the
  reported *r*/*R*<sup>2</sup>/*p* come from a single `lm()` fit across all
  taxa, and `geom_smooth()` is explicitly given `aes(group = 1)` so a
  single overall trend line (with one shared standard-error ribbon) is
  always drawn, rather than ggplot2 implicitly fitting one line per
  colour/shape level (which happens by default whenever colour/shape are
  mapped at the top level of `ggplot()` and inherited by `geom_smooth()`).
- **Statistic annotation placement.** The *r*/*R*<sup>2</sup>/*p* label is
  placed in the headroom created above the data (`y = Inf` with expanded
  upper plot margin, centred horizontally), so it never collides with
  points, error bars or the regression line regardless of whether the
  relationship is positive or negative.
- **Robust plotmath labels.** All parsed (`parse = TRUE`) plot annotations
  are built via `safe_plotmath_label()`, which validates the constructed
  expression with `parse()` before it reaches `ggplot2::annotate()`, and
  falls back to a plain-text label if a value is ever `NA`/non-finite,
  rather than letting a malformed label crash rendering.
- **Axis truncation.** Where one or two taxa have a disproportionately wide
  credible interval (e.g. `rate_genes1_median` for *H. influenzae* and
  *Neisseria* spp.), the axis is visually truncated with
  `coord_cartesian(xlim =, ylim =)` to keep the bulk of the points legible;
  this only crops the *view* — the regression line and reported statistics
  are always fitted to the complete, untruncated data, and this is noted
  explicitly in the affected figure's caption.
- **Legend layout.** Guide titles are placed *above* their (horizontal)
  row of keys (`guide_legend(title.position = "top")`) rather than to the
  left, which otherwise could make the combined legend wider than the
  figure.

## Data provenance and merge logic

- **Kwun & Croucher** data: taxon-level homologous recombination
  probabilities (from MLST-based generalised linear mixed model fits) and
  prophage density estimates (Poisson rate per isolate genome, from
  geNomad), together with transformability and ecological niche
  annotations, as reported in their Figure 2 and Table S2/S5, provided here
  as the raw per-taxon output of their own analysis pipeline.
- **Horsfield *et al.*** data: per-species PopPUNK-mod posterior medians for
  basal gene turnover rate, core mutation rate and proportion of fast
  genes (`Supplementary File 1`); pangenome size/composition statistics
  from WTBCluster; and negative-exponential core-vs-accessory distance
  curve fit parameters (Plateau, Rate ratio, Y-intercept) for a 22-species
  subset (`Supplementary File 8`).
- All merging (genus-level "spp." taxa matched to the mean of all
  corresponding species; *Enterococcus faecium* matched to
  `Enterococcus_B faecium`; taxa present in only one dataset excluded) is
  performed by `build_analysis_datasets()` in `R/data_prep.R`, executed at
  knit time — see `data/README.md` for the complete matching logic and a
  full account of which taxa are included/excluded at each merge step.
