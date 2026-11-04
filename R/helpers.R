# ---------------------------------------------------------------------------
# Helper functions, palettes and ggplot theming for the pangenome/HR/prophage
# conflict comparison analysis.
#
# Styling conventions (palettes, plotmath display labels, violin+jitter
# layout, label_parsed facets, geom_label_repel annotation of species) are
# deliberately kept consistent with
# https://github.com/nickjcroucher/recombination-conflict
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(scales)
  library(broom)
})

# ---- Palettes (as used in recombination-conflict) --------------------------

niche_palette <- c(
  "Environmental"           = 8,
  "Gastrointestinal tract"  = 16,
  "Respiratory tract"       = 1,
  "Multiple"                = 5
)

transformability_palette <- c(
  "Transformable"     = "cornflowerblue",
  "Non-transformable" = "purple"
)

# Shared minimal theme, matching the look of the source repository
theme_conflict <- function(...) {
  theme_bw(...) +
    theme(
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      legend.position   = "bottom"
    )
}

# ---- Lookup table of known uncertainty (95% CI) columns --------------------
#
# Several variables used across these comparisons have genuine, data-derived
# uncertainty intervals available in the merged CSVs (Poisson MLE CIs for
# lambda_hat; asymptotic CIs from the recombination-probability GLMM;
# Bayesian credible intervals from Horsfield et al.'s PopPUNK-mod fits).
# Where a lwr/upr pair is registered here, plot_correlation() will
# automatically draw the corresponding horizontal/vertical error bars.
# Variables not listed here (e.g. raw pangenome gene counts, genome size, and
# the negative-exponential Plateau/Rate-ratio/Y-intercept point estimates)
# have no reported uncertainty in the source data, so no bars are fabricated
# for them; this is flagged explicitly in the figure caption instead.

ci_lookup <- list(
  "lambda_hat"                  = c(lwr = "lambda_lwr",            upr = "lambda_upr"),
  "Recombination probability"   = c(lwr = "asymp.LCL",             upr = "asymp.UCL"),
  "rate_genes1_median"          = c(lwr = "rate_genes1_lower_CI",  upr = "rate_genes1_upper_CI"),
  "core_mu_median"               = c(lwr = "core_mu_lower_CI",      upr = "core_mu_upper_CI"),
  "prop_genes2_median"           = c(lwr = "prop_genes2_lower_CI",  upr = "prop_genes2_upper_CI")
)

#' Report, as a one-line markdown string, whether a variable has a registered
#' 95% CI, for use in figure captions.
ci_caption_note <- function(var_name) {
  if (var_name %in% names(ci_lookup)) {
    paste0("*", var_name, "* error bars show 95% CIs.")
  } else {
    paste0("*", var_name, "* has no reported uncertainty interval in the ",
          "source data, so no error bars are shown for this variable.")
  }
}

#' Safely coerce a single numeric value to a short, always-valid, non-empty
#' numeral string for embedding in a plotmath expression. Guards against
#' NULL/zero-length/NA/NaN/Inf values (e.g. from a degenerate model fit)
#' ever producing a malformed or truncated label string that would crash
#' `parse()` at render time.
safe_num <- function(x, digits = 2) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
  formatC(signif(x, digits), format = "g", digits = digits)
}

#' Build a plotmath annotation label from a template and named numeric
#' values, validating with parse() before use so that any unforeseen
#' formatting issue degrades gracefully to a plain-text (unparsed) label
#' rather than crashing ggplot's rendering step.
safe_plotmath_label <- function(template, ..., digits = 2) {
  label <- do.call(sprintf, c(list(template), lapply(list(...), safe_num, digits = digits)))
  parse_ok <- tryCatch({ parse(text = label); TRUE }, error = function(e) FALSE)
  list(label = label, parse = parse_ok)
}

# ---- Generic scatterplot with regression line + R^2/p annotation -----------
#
# Mirrors plot_rec_vs_element_model() from recombination-conflict, including
# geom_errorbar()/geom_errorbarh() for genuine 95% CIs/credible intervals
# where available (see ci_lookup above), but fits a simple linear model on
# (optionally log10-transformed) axes rather than a beta/GLMM, since here we
# are relating two continuous, taxon-level summary statistics rather than
# raw per-isolate counts.

plot_correlation <- function(df,
                             x_var, y_var,
                             x_lab = x_var, y_lab = y_var,
                             log_x = FALSE, log_y = FALSE,
                             colour_var = "Transformability",
                             shape_var  = "Niche",
                             label_var  = "display",
                             x_lwr_var = NULL, x_upr_var = NULL,
                             y_lwr_var = NULL, y_upr_var = NULL) {

  # Auto-detect CI columns from the lookup table unless explicitly overridden
  if (is.null(x_lwr_var) && x_var %in% names(ci_lookup)) x_lwr_var <- ci_lookup[[x_var]][["lwr"]]
  if (is.null(x_upr_var) && x_var %in% names(ci_lookup)) x_upr_var <- ci_lookup[[x_var]][["upr"]]
  if (is.null(y_lwr_var) && y_var %in% names(ci_lookup)) y_lwr_var <- ci_lookup[[y_var]][["lwr"]]
  if (is.null(y_upr_var) && y_var %in% names(ci_lookup)) y_upr_var <- ci_lookup[[y_var]][["upr"]]

  has_x_ci <- !is.null(x_lwr_var) && !is.null(x_upr_var) &&
    x_lwr_var %in% names(df) && x_upr_var %in% names(df)
  has_y_ci <- !is.null(y_lwr_var) && !is.null(y_upr_var) &&
    y_lwr_var %in% names(df) && y_upr_var %in% names(df)

  plot_df <- df %>%
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))

  x_fit <- if (log_x) log10(plot_df[[x_var]]) else plot_df[[x_var]]
  y_fit <- if (log_y) log10(plot_df[[y_var]]) else plot_df[[y_var]]

  fit <- stats::lm(y_fit ~ x_fit)
  fit_summary <- broom::glance(fit)
  r2     <- fit_summary$r.squared
  p_val  <- fit_summary$p.value
  slope  <- stats::coef(fit)[["x_fit"]]
  r_val  <- if (is.null(slope) || is.na(slope) || is.na(r2)) NA_real_ else sign(slope) * sqrt(r2)

  stat_label_res <- safe_plotmath_label(
    "italic(r) == %s*'; '*italic(R)^2 == %s*'; '*italic(p) == %s",
    r_val, r2, p_val
  )
  stat_label <- stat_label_res$label
  stat_label_parse <- stat_label_res$parse

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      colour = .data[[colour_var]],
      shape  = .data[[shape_var]],
      label  = .data[[label_var]]
    )
  ) +
    ggplot2::geom_smooth(
      method = "lm", se = TRUE, colour = "black", fill = "grey70",
      linewidth = 0.8, alpha = 0.3
    )

  # Vertical error bars (uncertainty in y)
  if (has_y_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[y_lwr_var]], ymax = .data[[y_upr_var]]),
      width = 0, alpha = 0.5
    )
  }

  # Horizontal error bars (uncertainty in x)
  if (has_x_ci) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data[[x_lwr_var]], xmax = .data[[x_upr_var]]),
      height = 0, alpha = 0.5
    )
  }

  p <- p +
    ggplot2::geom_point(size = 2.5) +
    ggrepel::geom_label_repel(size = 3, force_pull = 0.5, parse = TRUE,
                              show.legend = FALSE, max.overlaps = 30) +
    ggplot2::scale_shape_manual(values = niche_palette) +
    ggplot2::scale_colour_manual(values = transformability_palette,
                                 aesthetics = c("colour", "fill")) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::annotate(
      "text",
      x = unname(stats::quantile(plot_df[[x_var]], 0.98, na.rm = TRUE)),
      y = unname(stats::quantile(plot_df[[y_var]], 0.98, na.rm = TRUE)),
      label = stat_label, parse = stat_label_parse, hjust = 1, size = 3.5
    ) +
    theme_conflict()

  if (log_x) p <- p + ggplot2::scale_x_log10()
  if (log_y) p <- p + ggplot2::scale_y_log10()

  caption_note <- paste(ci_caption_note(x_var), ci_caption_note(y_var))

  list(plot = p, model = fit, r = r_val, r2 = r2, p_value = p_val,
      data = plot_df, has_x_ci = has_x_ci, has_y_ci = has_y_ci,
      caption_note = caption_note)
}

# ---- Violin + jitter comparison across a binary/categorical grouping -------
#
# Mirrors plot_violin_kw() from recombination-conflict: a Wilcoxon rank-sum
# test annotated in the top-right corner, points coloured by Transformability
# and shaped by Niche, median line drawn through the violin.

plot_violin_comparison <- function(df, x_var, y_var,
                                   x_lab = x_var, y_lab = y_var,
                                   colour_var = "Transformability",
                                   shape_var  = "Niche") {

  plot_df <- df %>%
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))

  x_vec <- plot_df[[x_var]]
  y_vec <- plot_df[[y_var]]

  test_res <- stats::wilcox.test(y_vec ~ x_vec)
  p_label_res <- safe_plotmath_label(
    "Wilcoxon~rank~sum~test~italic(p) == %s", test_res$p.value, digits = 3
  )
  p_label <- p_label_res$label
  p_label_parse <- p_label_res$parse

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_violin(draw_quantiles = 0.5, alpha = 0.3, colour = "grey80") +
    ggplot2::geom_jitter(
      ggplot2::aes(colour = .data[[colour_var]], shape = .data[[shape_var]]),
      width = 0.15, height = 0, size = 2.5, alpha = 0.85
    ) +
    ggplot2::scale_shape_manual(values = niche_palette) +
    ggplot2::scale_colour_manual(values = transformability_palette) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf, label = p_label,
      hjust = 1.05, vjust = 1.25, size = 3.5, parse = p_label_parse
    ) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    theme_conflict() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  list(plot = p, test = test_res, data = plot_df)
}

# ---- Combine individual panels + shared legend, as in the summary figure ---

combine_with_legend <- function(plotlist, labels = "AUTO", nrow = 2, ncol = 2,
                                legend_source = 1) {
  legend <- cowplot::get_plot_component(
    plotlist[[legend_source]] + ggplot2::theme(legend.position = "bottom"),
    "guide-box-bottom", return_all = TRUE
  )

  panels <- cowplot::plot_grid(
    plotlist = lapply(plotlist, function(p) p + ggplot2::theme(legend.position = "none")),
    nrow = nrow, ncol = ncol, labels = labels
  )

  cowplot::plot_grid(panels, legend, nrow = 2, rel_heights = c(0.92, 0.08))
}
