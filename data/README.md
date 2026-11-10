# Data dictionary

This repository stores only the **three raw, unmodified source files**
below in `data/raw/`. All taxon matching and merging is performed live by
`R/data_prep.R::build_analysis_datasets()` when the main R Markdown
document is knitted — no pre-merged CSV is committed to the repository, so
the full provenance of every number in the analysis is traceable back to
these three files.

## `raw/recombination_conflict_comparison_plot.csv`

Kwun & Croucher's own per-taxon output (n = 23 taxa), produced directly by
their
[recombination-conflict](https://github.com/nickjcroucher/recombination-conflict)
pipeline (the `prophage_rec_df` object, filtered to the geNomad
Caudoviricetes detector/element combination).

| Column | Description |
|---|---|
| `species`, `display`, `display_name`, `Species` | Taxon name (full binomial / plotmath-formatted / short code) |
| `detector`, `element_type` | `geNomad`, `Caudoviricetes` |
| `n_isolates`, `k_total` | Number of isolate genomes assayed and total prophage count |
| `lambda_hat`, `lambda_lwr`, `lambda_upr` | Poisson MLE (+ 95% CI) of the expected number of prophage per isolate genome |
| `mean_empirical`, `var_empirical`, `aic`, `bic`, `fit_converged` | Poisson model fit diagnostics |
| `Recombination probability`, `SE`, `df`, `asymp.LCL`, `asymp.UCL` | Probability of diversification through homologous recombination (vs point mutation) at a locus, and its asymptotic 95% CI, from the MLST-based beta-binomial GLMM |
| `Transformability` | `Transformable` / `Non-transformable`, from laboratory evidence |
| `Niche` | Primary ecological niche |

## `raw/supplementary_file1_408species.csv`

Horsfield *et al.*'s **Supplementary File 1**: one row per bacterial
species (n = 408), giving taxonomy (GTDB), PopPUNK-mod pangenome-dynamics
parameter estimates, and pangenome composition statistics computed with
WTBCluster. Key columns used here:

| Column | Description |
|---|---|
| `Genus`, `Species` | GTDB taxonomy (note some genera are suffixed, e.g. `Enterococcus_B`, to disambiguate polyphyletic groups) |
| `rate_genes1_median` (+ `_lower_CI`/`_upper_CI`) | PopPUNK-mod posterior median basal gene turnover rate |
| `core_mu_median` (+ CIs) | PopPUNK-mod posterior median core mutation rate |
| `prop_genes2_median` (+ CIs) | PopPUNK-mod posterior median proportion of "fast"-exchanged accessory genes |
| `generalism_score`, `habitat_count` | Ecological generalism score (metaTraits) and number of habitats occupied |
| `pangenome_fluidity`, `pan_genes`, `core_genes`, `intermediate_genes`, `rare_genes`, `avg_gene_freq` | Pangenome composition statistics (WTBCluster) |
| `percentage_core_80/90/100`, `accessory_80/90/100`, `core_80/90/100` | Core/accessory genome size at varying frequency thresholds |
| `pan_size (Dewar et al)`, `Ne`, `genome_size`, `Number_of_Strain` | Additional pangenome and effective population size estimates |

Many ecological/habitat columns (`Host_or_free`, `Primary_environment`,
`Mammal`, `Soil`, ... `Total_envs`) are populated only for the subset of
species with a `generalism_score`; they are not used in this comparison
analysis.

## `raw/supplementary_file8_negexp_fits.csv`

Horsfield *et al.*'s **Supplementary File 8**: negative-exponential
core-vs-accessory pairwise distance curve fits (their Equation 1,
*y = i − (i − j)e^(−kx)*), in long format, for 22 species/datasets
(abbreviated binomials, e.g. `E. coli`; *S. pneumoniae* has separate
`(Global)` and `(Local)` fits).

| Column | Description |
|---|---|
| `Species` | Abbreviated binomial / dataset label |
| `Coefficient` | `Y-intercept (j)`, `Rate ratio (k)`, or `Plateau (i)` |
| `Estimate`, `Upper 95% CI`, `Lower CI`, `p-value` | Fitted value, 95% CI bounds and significance |

## Taxon matching logic (`R/data_prep.R`)

**`main_df` (n = 23):** every taxon in
`recombination_conflict_comparison_plot.csv` is matched to Supplementary
File 1 as follows:

- **Direct match (16 taxa):** `species` matches a single Supplementary
  File 1 `Species` value exactly (e.g. `Escherichia coli`).
- **Genus-level mean (6 taxa):** *Campylobacter* spp., *Citrobacter*
  spp., *Enterobacter* spp., *Neisseria* spp., *Proteus* spp. and
  *Serratia* spp. are matched to the unweighted mean of all Supplementary
  File 1 species sharing that `Genus`.
- **Renamed match (1 taxon):** *Enterococcus faecium* is matched to the
  GTDB-style label `Enterococcus_B faecium` in Supplementary File 1.

The number of Supplementary File 1 rows contributing to each taxon's
values is recorded in `n_species_matched`.

**`negexp_df` (n = 15, subset of `main_df`):** Supplementary File 8 uses a
different abbreviated-binomial naming convention and covers a
non-identical set of species, so only 15 of the 23 taxa have a match:

- **Direct match (12 taxa):** *A. baumanii* → *Acinetobacter baumannii*,
  *E. coli* → *Escherichia coli*, *E. faecalis* → *Enterococcus
  faecalis*, *H. influenzae* → *Haemophilus influenzae*, *H. pylori* →
  *Helicobacter pylori*, *K. pneumoniae* → *Klebsiella pneumoniae*,
  *M. catarrhalis* → *Moraxella catarrhalis*, *P. aeruginosa* →
  *Pseudomonas aeruginosa*, *S. enterica* → *Salmonella enterica*,
  *S. aureus* → *Staphylococcus aureus*, *S. agalactiae* →
  *Streptococcus agalactiae*, *S. pyogenes* → *Streptococcus pyogenes*.
- **Averaged match (1 taxon):** *Streptococcus pneumoniae* = mean of the
  `S. pneumoniae (Global)` and `S. pneumoniae (Local)` fits.
- **Genus-level match (2 taxa):** *Neisseria* spp. = mean of
  *N. gonorrhoeae* and *N. meningitidis*; *Campylobacter* spp. =
  *C. jejuni* only (the sole *Campylobacter* species fitted in
  Supplementary File 8).
- **No match (8 taxa, excluded from `negexp_df`):** *Bacillus subtilis*,
  *Citrobacter* spp., *Enterobacter* spp., *Enterococcus faecium*,
  *Morganella morganii*, *Proteus* spp., *Serratia* spp., *Vibrio
  cholerae*.

`Y-intercept x genome size` is computed on the fly in `negexp_df` as
`` `Y-intercept (j)` * genome_size ``; *Moraxella catarrhalis* lacks a
`genome_size` value in Supplementary File 1 and so resolves to `NA` for
this derived quantity (n = 14 for comparisons using it).

### Note on duplicated species rows in Supplementary File 1

Five species (`Bacillus_A cereus`, `Bacillus_A thuringiensis_S`,
`Enterobacter asburiae_B`, `Haemophilus influenzae`, `Serratia
sarumanii`) each have two rows in Supplementary File 1, with distinct
`tip` (representative genome) identifiers but an identical `Species`
label — reflecting two independently fitted representative genomes for
the same GTDB species. `match_horsfield_to_taxa()` treats this exactly
like a genus-level match: it averages across *all* rows sharing that
`Species` label, so *Haemophilus influenzae*'s values in `main_df`
(`n_species_matched = 2`) are the mean of its two Supplementary File 1
rows, and the duplicate *Serratia sarumanii* / *Enterobacter asburiae_B*
rows are both included (each counted once) within the *Serratia* spp. /
*Enterobacter* spp. genus-level means.

## Sample size calculations

The counts of isolates in `supplementary_file12_negExp_isolate_counts` and `supplementary_file12_isolate_counts`
were calculated from the original Supplementary File 12 of Horsfield *et al*. and stored in summarised
form to keep the repository size small. The processing used the workflow:

```
supp12_df %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n_isolates = dplyr::n_distinct(Accession))
```