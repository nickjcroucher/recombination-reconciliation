# ---------------------------------------------------------------------------
# Data preparation: build the two analysis datasets directly from the three
# raw, unmodified source files in data/raw/:
#
#   1. supplementary_file1_408species.csv
#        Horsfield et al.'s Supplementary File 1 -- PopPUNK-mod pangenome
#        dynamics parameters and pangenome composition statistics for 408
#        bacterial species (one row per species).
#
#   2. supplementary_file8_negexp_fits.csv
#        Horsfield et al.'s Supplementary File 8 -- negative-exponential
#        core-vs-accessory pairwise distance curve fits (Equation 1:
#        y = i - (i - j)e^(-kx)) for 22 species/datasets, in long format
#        (one row per Species x Coefficient).
#
#   3. recombination_conflict_comparison_plot.csv
#        Kwun & Croucher's per-taxon output (n = 23 taxa) from their
#        https://github.com/nickjcroucher/recombination-conflict pipeline:
#        geNomad Caudoviricetes prophage density (lambda_hat, Poisson MLE +
#        95% CI) and the probability of diversification through homologous
#        recombination (Recombination probability, from the MLST-based
#        beta-binomial GLMM + asymptotic CI), together with Transformability
#        and Niche annotations.
#
# No pre-merged CSV is stored in the repository: build_analysis_datasets()
# performs the full taxon-matching/aggregation logic below every time the
# Rmd is knitted, so the merge is transparent and reproducible from the
# original supplementary files alone.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

# Numeric Horsfield et al. columns (Supplementary File 1) carried through
# into the merged dataset.
horsfield_vars <- c(
  "rate_genes1_median", "rate_genes1_lower_CI", "rate_genes1_upper_CI",
  "core_mu_median",     "core_mu_lower_CI",     "core_mu_upper_CI",
  "prop_genes2_median", "prop_genes2_lower_CI", "prop_genes2_upper_CI",
  "generalism_score",   "habitat_count",        "pangenome_fluidity",
  "percentage_core_100", "percentage_core_90",  "percentage_core_80",
  "percentage_accessory_20", "percentage_accessory_10",
  "pan_size (Dewar et al)", "accessory_90", "accessory_80",
  "Ne", "genome_size", "Number_of_Strain",
  "pan_genes", "core_genes", "intermediate_genes", "rare_genes",
  "avg_gene_freq"
)

# ---- Taxon matching between Kwun & Croucher's 23 taxa and Horsfield et
# al.'s Supplementary File 1 ------------------------------------------------
#
# Most of the 23 taxa match a single species' row directly by binomial name.
# Genus-level ("spp.") taxa are matched by taking the unweighted mean across
# all matching species present in Supplementary File 1 for that genus.
# Enterococcus faecium is recorded under the GTDB-style label
# "Enterococcus_B faecium" in Supplementary File 1.

genus_mean_taxa <- c(
  "Campylobacter spp" = "Campylobacter",
  "Citrobacter spp"   = "Citrobacter",
  "Enterobacter spp"  = "Enterobacter",
  "Neisseria spp"     = "Neisseria",
  "Proteus spp"       = "Proteus",
  "Serratia spp"      = "Serratia"
)

species_rename <- c(
  "Enterococcus faecium" = "Enterococcus_B faecium"
)

#' Match each of Kwun & Croucher's 23 taxa to Horsfield et al.'s
#' Supplementary File 1, averaging across species for genus-level taxa.
match_horsfield_to_taxa <- function(taxa, file1_df) {

  purrr::map_dfr(taxa, function(taxon) {

    if (taxon %in% names(genus_mean_taxa)) {
      # Genus-level ("spp.") taxon: mean across all matching species
      genus <- genus_mean_taxa[[taxon]]
      matched <- file1_df %>% dplyr::filter(Genus == genus)
      n_matched <- nrow(matched)

      summary_row <- matched %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(horsfield_vars),
                                       ~ mean(.x, na.rm = TRUE)))

    } else {
      lookup_species <- if (taxon %in% names(species_rename)) {
        unname(species_rename[[taxon]])
      } else {
        taxon
      }
      matched <- file1_df %>% dplyr::filter(Species == lookup_species)
      n_matched <- nrow(matched)

      summary_row <- matched %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(horsfield_vars),
                                       ~ mean(.x, na.rm = TRUE)))
    }

    summary_row %>%
      dplyr::mutate(species = taxon, n_species_matched = n_matched,
                    .before = 1)
  })
}

# ---- Taxon matching between Kwun & Croucher's 23 taxa and Horsfield et
# al.'s Supplementary File 8 (negative-exponential curve fits) -------------
#
# Supplementary File 8 uses abbreviated binomials (e.g. "E. coli") for a
# different set of 22 species/datasets. Only 15 of the 23 taxa have a
# match: 12 direct one-to-one matches, S. pneumoniae (averaged across the
# "Global" and "Local" fits), and two genus-level matches (Campylobacter
# spp -> C. jejuni only; Neisseria spp -> mean of N. gonorrhoeae and
# N. meningitidis).

negexp_direct_map <- c(
  "A. baumanii"    = "Acinetobacter baumannii",
  "E. coli"        = "Escherichia coli",
  "E. faecalis"    = "Enterococcus faecalis",
  "H. influenzae"  = "Haemophilus influenzae",
  "H. pylori"      = "Helicobacter pylori",
  "K. pneumoniae"  = "Klebsiella pneumoniae",
  "M. catarrhalis" = "Moraxella catarrhalis",
  "P. aeruginosa"  = "Pseudomonas aeruginosa",
  "S. enterica"    = "Salmonella enterica",
  "S. aureus"      = "Staphylococcus aureus",
  "S. agalactiae"  = "Streptococcus agalactiae",
  "S. pyogenes"    = "Streptococcus pyogenes"
)

negexp_averaged_map <- list(
  "Streptococcus pneumoniae" = c("S. pneumoniae (Global)", "S. pneumoniae (Local)"),
  "Neisseria spp"            = c("N. gonorrhoeae", "N. meningitidis"),
  "Campylobacter spp"        = c("C. jejuni")
)

#' Pivot Supplementary File 8 to one row per File 8 "species" label, with
#' Plateau (i), Rate ratio (k) and Y-intercept (j) as columns.
widen_negexp <- function(file8_df) {
  file8_df %>%
    dplyr::select(Species, Coefficient, Estimate) %>%
    tidyr::pivot_wider(names_from = Coefficient, values_from = Estimate)
}

#' Match Kwun & Croucher taxa to Horsfield et al.'s negative-exponential
#' curve fits (Supplementary File 8), returning only the taxa with a match.
match_negexp_to_taxa <- function(negexp_wide) {

  direct_df <- tibble::tibble(
    file8_species = names(negexp_direct_map),
    species        = unname(negexp_direct_map)
  ) %>%
    dplyr::inner_join(negexp_wide, by = c("file8_species" = "Species")) %>%
    dplyr::select(-file8_species)

  averaged_df <- purrr::imap_dfr(negexp_averaged_map, function(labels, taxon) {
    negexp_wide %>%
      dplyr::filter(Species %in% labels) %>%
      dplyr::summarise(
        `Plateau (i)`    = mean(`Plateau (i)`, na.rm = TRUE),
        `Rate ratio (k)` = mean(`Rate ratio (k)`, na.rm = TRUE),
        `Y-intercept (j)`= mean(`Y-intercept (j)`, na.rm = TRUE)
      ) %>%
      dplyr::mutate(species = taxon, .before = 1)
  })

  dplyr::bind_rows(direct_df, averaged_df)
}

# ---- Top-level entry point -------------------------------------------------

#' Build the two taxon-level analysis datasets directly from the three raw
#' source files.
#'
#' @return a list with elements `main_df` (n = 23 taxa: Kwun & Croucher's HR
#'   probability and prophage density joined to Horsfield et al.'s
#'   pangenome-dynamics parameters) and `negexp_df` (the n = 15 subset of
#'   `main_df` additionally joined to the negative-exponential curve fit
#'   parameters, with `Y-intercept x genome size` pre-computed).
build_analysis_datasets <- function(
    file1_path  = "data/raw/supplementary_file1_408species.csv",
    file8_path  = "data/raw/supplementary_file8_negexp_fits.csv",
    plot_path   = "data/raw/recombination_conflict_comparison_plot.csv"
) {

  rec_df <- readr::read_csv(plot_path, show_col_types = FALSE)

  file1_df <- readr::read_csv(file1_path, show_col_types = FALSE)

  file8_df <- readr::read_csv(file8_path, show_col_types = FALSE)

  # --- Merge 1: Kwun & Croucher taxa <- Horsfield et al. pangenome params
  horsfield_df <- match_horsfield_to_taxa(rec_df$species, file1_df)

  main_df <- rec_df %>%
    dplyr::left_join(horsfield_df, by = "species") %>%
    dplyr::mutate(
      Transformability = factor(Transformability,
                                levels = c("Transformable", "Non-transformable")),
      Niche = factor(Niche)
    )

  # --- Merge 2: subset of taxa <- Horsfield et al. negative-exponential fits
  negexp_wide <- widen_negexp(file8_df)
  negexp_match_df <- match_negexp_to_taxa(negexp_wide)

  negexp_df <- main_df %>%
    dplyr::inner_join(negexp_match_df, by = "species") %>%
    dplyr::mutate(`Y-intercept x genome size` = `Y-intercept (j)` * genome_size)

  list(main_df = main_df, negexp_df = negexp_df)
}
