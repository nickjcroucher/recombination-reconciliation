# ---------------------------------------------------------------------------
# R/data_prep.R
#
# Builds every analysis dataset used in pangenome_recombination_conflict.Rmd
# directly from the three raw input files in data/raw/. No merged/derived
# CSVs are stored in the repository - everything is recomputed on each knit.
#
# Inputs (data/raw/):
#   - supplementary_file1_408species.csv        Horsfield et al., Supplementary File 1
#                                                (408-species PopPUNK-mod parameter estimates
#                                                 and pangenome statistics)
#   - supplementary_file8_negexp_fits.csv        Horsfield et al., Supplementary File 8
#                                                (negative-exponential curve fits: Y-intercept,
#                                                 Plateau, Rate ratio, for 22 species + 1 extra
#                                                 S. pneumoniae local-sample fit)
#   - recombination_conflict_comparison_plot.csv Kwun & Croucher, Figure 2 source data
#                                                (23 taxa: HR probability from MLST, and
#                                                 geNomad-detected intact Caudoviricetes
#                                                 prophage counts per isolate, fitted as a
#                                                 Poisson rate lambda_hat)
#
# IMPORTANT SCOPE NOTE: "prophage" throughout this analysis refers specifically to
# geNomad-detected, hallmark-gene-rich Caudoviricetes provirus regions - i.e. large,
# likely-intact temperate phage capable of causing lytic induction and cell death.
# Satellite phage/phage-related chromosomal islands (e.g. PRCIs), inoviruses, and
# plasmids were explicitly excluded from Kwun & Croucher's cross-taxon genomic
# comparison (Fig. 2A) and are therefore also excluded here.
# ---------------------------------------------------------------------------

#' Build a lookup from a Species binomial/"spp" label (as used in the
#' recombination-conflict comparison dataset) to a species or genus label as
#' used in Supplementary File 1 (Horsfield et al.)
genus_level_taxa <- c(
  "Campylobacter spp", "Citrobacter spp", "Enterobacter spp",
  "Neisseria spp", "Proteus spp", "Serratia spp"
)

# Enterococcus faecium is recorded under the GTDB-style genus label
# "Enterococcus_B" in Supplementary File 1
species_rename <- list(
  "Enterococcus faecium" = "Enterococcus_B faecium"
)

#' Match each of the 23 taxa in the HR/prophage comparison dataset to
#' Supplementary File 1, averaging across member species for genus-level
#' ("spp") taxa, and renaming where GTDB taxonomy differs from the binomial
#' used by Kwun & Croucher.
match_horsfield_to_taxa <- function(taxon_vec, file1_df,
                                    return_matched_species = FALSE,
                                    return_sum = FALSE) {
  
  numeric_cols <- names(file1_df)[vapply(file1_df, is.numeric, logical(1))]
  
  purrr::map_dfr(unique(taxon_vec), function(taxon) {
    
    is_genus_level <- taxon %in% genus_level_taxa
    
    if (is_genus_level) {
      genus <- sub(" spp$", "", taxon)
      matched_rows <- file1_df %>%
        dplyr::filter(Genus == genus)
    } else {
      if (taxon %in% names(species_rename)) {
        lookup_species <- unname(species_rename[[taxon]])
      } else {
        lookup_species <- taxon
      }
      
      matched_rows <- file1_df %>%
        dplyr::filter(Species == lookup_species)
    }
    
    if (nrow(matched_rows) == 0) {
      result <- tibble::tibble(
        species = taxon,
        n_matched = 0L
      )
      
      if (return_matched_species) {
        result$matched_species <- NA_character_
      }
      
      return(result)
    }
    
    summarised <- NULL
    if (return_sum) {
      summarised <- matched_rows %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(numeric_cols),
            ~ sum(.x, na.rm = TRUE)
          )
        )
    } else {
      summarised <- matched_rows %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(numeric_cols),
            ~ mean(.x, na.rm = TRUE)
          )
        )
    }

    result <- dplyr::bind_cols(
      tibble::tibble(
        species = taxon,
        n_matched = nrow(matched_rows)
      ),
      summarised
    )
    
    if (return_matched_species) {
      result <- tibble::add_column(
        result,
        matched_species = paste(unique(matched_rows$Species), collapse = "; "),
        .after = "species"
      )
    }
    
    result
  })
}

#' Match the 23-taxon comparison dataset to Supplementary File 8
#' (negative-exponential fits). Only a subset of taxa have a corresponding
#' fit; "spp" taxa are matched by averaging across any member species
#' present, and S. pneumoniae is averaged across its Global/Local fits.
negexp_species_map <- list(
  "Acinetobacter baumannii"   = c("A. baumanii"),
  "Bacillus subtilis"         = c(),  # no File 8 entry
  "Campylobacter spp"         = c("C. jejuni", "Campylobacter jejuni"),
  "Citrobacter spp"           = c(),
  "Enterobacter spp"          = c(),
  "Enterococcus faecalis"     = c("E. faecalis"),
  "Enterococcus faecium"      = c(),
  "Escherichia coli"          = c("E. coli"),
  "Haemophilus influenzae"    = c("H. influenzae"),
  "Helicobacter pylori"       = c("H. pylori"),
  "Klebsiella pneumoniae"     = c("K. pneumoniae"),
  "Moraxella catarrhalis"     = c("M. catarrhalis", "Morexalla catarrhalis"),
  "Morganella morganii"       = c(),
  "Neisseria spp"             = c("N. gonorrhoeae", "N. meningitidis", "Neisseria gonorrhoeae", "Neisseria meningitidis"),
  "Proteus spp"               = c(),
  "Pseudomonas aeruginosa"    = c("P. aeruginosa"),
  "Salmonella enterica"       = c("S. enterica"),
  "Serratia spp"              = c(),
  "Staphylococcus aureus"     = c("S. aureus"),
  "Streptococcus agalactiae"  = c("S. agalactiae"),
  "Streptococcus pneumoniae"  = c("S. pneumoniae (Global)"),
  "Streptococcus pyogenes"    = c("S. pyogenes"),
  "Vibrio cholerae"           = c()
)

match_negexp_to_taxa <- function(taxon_vec, file8_df) {

  file8_wide <- file8_df %>%
    dplyr::filter(Species_short != "S. pneumoniae (Local)") %>% 
    dplyr::select(Species, Parameter_name, Median, Upper_CI, Lower_CI) %>%
    tidyr::pivot_wider(
      names_from  = Parameter_name,
      values_from = c(Median, Upper_CI, Lower_CI),
      names_glue  = "{Parameter_name}_{.value}"
    )

  numeric_cols <- names(file8_wide)[vapply(file8_wide, is.numeric, logical(1))]

  purrr::map_dfr(unique(taxon_vec), function(taxon) {

    member_species <- negexp_species_map[[taxon]]
    if (is.null(member_species) || length(member_species) == 0) {
      return(NULL)
    }

    matched_rows <- file8_wide %>% dplyr::filter(Species %in% member_species | Species %in% c(taxon))
    if (nrow(matched_rows) == 0) return(NULL)

    summarised <- matched_rows %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)))

    dplyr::bind_cols(tibble::tibble(species = taxon), summarised)
  })
}

get_sample_sizes <- function(df1,df2,species_list) {
  
  df1 %>% 
    dplyr::rename(Kwun_n_isolates = any_of("n_isolates")) %>% 
    dplyr::left_join(match_horsfield_to_taxa(species_list,
                                             df2 %>%
                                               tidyr::separate_wider_delim(Species,
                                                                           names=c("Genus","specific_name"),
                                                                           delim=" ",
                                                                           cols_remove = FALSE,
                                                                           too_many = "drop") %>%
                                               dplyr::select(Species,Genus,n_isolates) %>%
                                               dplyr::distinct(),
                                             return_sum = TRUE
    ),
    by = "species"
    ) %>% 
    dplyr::rename(Horsfield_n_isolates = any_of("n_isolates"))
  
}

standardise_species <- function(x) {
  
  x <- stringr::str_replace(
    x,
    "^([A-Za-z]+)_([A-Z]) ",
    function(m) {
      parts <- stringr::str_match(
        m, "^([A-Za-z]+)_([A-Z]) "
      )
      paste0(parts[, 2], " ", tolower(parts[, 3]), " ")
    }
  )
  
  x <- stringr::str_remove(x, "_[A-Z]+$")
  x <- stringr::str_remove(x, " [a-z]$")
  
  x
}

#' Master function: builds both analysis datasets from the three raw files.
#'
#' Returns a list with:
#'   full_df   - the 23-taxon HR/prophage dataset merged with Supplementary
#'               File 1 pangenome/PopPUNK-mod parameters
#'   negexp_df - the subset (n=15) of full_df additionally merged with
#'               Supplementary File 8 negative-exponential fit parameters,
#'               plus a genome-size-scaled Y-intercept column
#'   raw_atb  - the model fits to the 408 species dataset
build_analysis_datasets <- function(raw_dir = "data/raw") {
  
  # Parse data files
  rec_df   <- readr::read_csv(file.path(raw_dir, "recombination_conflict_comparison_plot.csv"),
                              show_col_types = FALSE)
  file8_df <- readr::read_csv(file.path(raw_dir, "supplementary_file8_negexp_fits.csv"),
                              show_col_types = FALSE) %>%
    dplyr::mutate(Species =
                    dplyr::if_else(
                      Species=="Morexalla catarrhalis",
                      "Moraxella catarrhalis",
                      Species)
    )
  file1_df <- readr::read_csv(file.path(raw_dir, "supplementary_file1_408species.csv"),
                              show_col_types = FALSE) %>%
    dplyr::mutate(Species =
                    dplyr::if_else(
                      Species=="Morexalla catarrhalis",
                      "Moraxella catarrhalis",
                      Species)
    )
  full_sample_size_df <-
    readr::read_csv("./data/raw/supplementary_file12_isolate_counts.csv")
  negexp_sample_size_df <-
    readr::read_csv("./data/raw/supplementary_file12_negExp_isolate_counts.csv")
  
  # Match taxa names
  horsfield_matched <-
    match_horsfield_to_taxa(rec_df$species, file1_df) %>%
    dplyr::left_join(
      match_horsfield_to_taxa(rec_df$species,
                              file8_df %>%
                                tidyr::separate_wider_delim(Species,
                                                            names=c("Genus","specific_name"),
                                                            delim=" ",
                                                            cols_remove = FALSE,
                                                            too_many = "drop") %>%
                                dplyr::mutate(Mean_genome_size = file8_df$Pangenome_size*file8_df$Avg_genome_size) %>%
                                dplyr::select(Species,Genus,Mean_genome_size) %>%
                                dplyr::distinct()
      ) %>% dplyr::select(-n_matched),
      by = "species"
    )
  
  full_df <- rec_df %>%
    dplyr::left_join(horsfield_matched,
                     by = "species")
  
  negexp_matched <- match_negexp_to_taxa(rec_df$species, file8_df)
  
  # Process rates
  negexp_df <- full_df %>%
    dplyr::inner_join(negexp_matched, by = "species") %>%
    dplyr::mutate(
      yint_times_genome       = `Proportion fast genes_Median`       * Mean_genome_size,
      yint_times_genome_lwr   = `Proportion fast genes_Lower_CI`       * Mean_genome_size,
      yint_times_genome_upr   = `Proportion fast genes_Upper_CI`   * Mean_genome_size
    )
  
  # Add sample sizes
  file1_numeric_cols <- names(file1_df)[vapply(file1_df, is.numeric, logical(1))]
  file1_df <-
    file1_df %>%
    dplyr::mutate(Species_std = standardise_species(Species)) %>% 
    dplyr::group_by(Species_std) %>% 
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(file1_numeric_cols),
        ~ sum(.x, na.rm = TRUE)
      )
    ) %>% 
    dplyr::group_by() %>% 
    dplyr::left_join(full_sample_size_df %>%
                       dplyr::mutate(Species_std = standardise_species(Species)) %>% 
                       dplyr::select(-Species),
                     by = "Species_std")
  full_df <- full_df %>% get_sample_sizes(full_sample_size_df,rec_df$species)
  negexp_df <- negexp_df %>% get_sample_sizes(negexp_sample_size_df,rec_df$species)
  
  # Return file datasets
  list(raw_atb = file1_df, full_df = full_df, negexp_df = negexp_df)
}