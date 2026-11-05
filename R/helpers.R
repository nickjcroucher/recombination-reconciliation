# ---------------------------------------------------------------------------
# R/helpers.R
#
# Shared palettes, theme, label-safety utilities, and plotting functions for
# pangenome_recombination_conflict.Rmd. Styled to match
# https://github.com/nickjcroucher/recombination-conflict
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(broom)
  library(scales)
})

# -----------------------------------------------------------------------
# Palettes (matching nickjcroucher/recombination-conflict conventions)
# -----------------------------------------------------------------------

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

theme_pangenome <- function() {
  theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      legend.box       = "vertical"
    )
}

# A single, shared horizontal-legend specification: titles sit ABOVE the
# keys (title.position = "top") and keys are laid out in one row, so the
# combined legend never becomes wider than the plot panel.
guides_pangenome <- function() {
  guides(
    colour = guide_legend(title.position = "top", nrow = 1, title.hjust = 0.5),
    fill   = guide_legend(title.position = "top", nrow = 1, title.hjust = 0.5),
    shape  = guide_legend(title.position = "top", nrow = 1, title.hjust = 0.5)
  )
}

# -----------------------------------------------------------------------
# Safe numeric formatting / plotmath label construction
#
# These exist because plotmath annotation strings built with naive
# paste0()/sprintf() calls fail in two distinct ways:
#   (1) if any input is NA/NaN/Inf/length-0, the resulting string is
#       malformed and parse() errors out only when ggplot actually tries
#       to draw the layer (i.e. far downstream of the real problem);
#   (2) chaining more than one relational operator (==) in a single
#       plotmath expression without parentheses is invalid R syntax, e.g.
#       `italic(r) == 0.85 * '; ' * italic(R)^2 == 0.72` does not parse,
#       because == is non-associative. Each "name == value" clause must be
#       individually parenthesised before being joined with an
#       associative operator (`*`).
# -----------------------------------------------------------------------

safe_num <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0) return("NA")
  x <- x[1]
  if (is.na(x) || is.nan(x) || is.infinite(x)) return("NA")
  format(signif(x, digits), trim = TRUE, scientific = abs(x) < 1e-4 || abs(x) >= 1e5)
}

#' Build a single parenthesised plotmath relation, e.g. "(italic(r) == 0.85)"
plotmath_relation <- function(name_expr, value, digits = 3) {
  paste0("(", name_expr, " == ", safe_num(value, digits), ")")
}

#' Combine multiple "name == value" relations into one valid, parse()-safe
#' plotmath expression by joining individually-parenthesised relations with
#' the (associative) multiplication operator, which plotmath renders as a
#' plain separator when interleaved with a quoted string.
plotmath_stat_relations <- function(relations, sep = "'; '") {
  paste(relations, collapse = paste0(" * ", sep, " * "))
}

#' Build a stats annotation label and verify it parses before returning it.
#' Falls back to an unparsed, informative text string (with parse = FALSE)
#' if validation fails for any reason, so a single bad value can never
#' crash the whole render.
safe_plotmath_label <- function(relations, sep = "'; '") {
  label <- plotmath_stat_relations(relations, sep = sep)
  parse_ok <- tryCatch({
    parse(text = label)
    TRUE
  }, error = function(e) FALSE)

  if (parse_ok) {
    list(label = label, parse = TRUE)
  } else {
    plain <- gsub("italic\\(([a-zA-Z0-9]+)\\)", "\\1", label)
    plain <- gsub("[()*]", " ", plain)
    plain <- gsub("\\^2", "2", plain)
    list(label = plain, parse = FALSE)
  }
}

# -----------------------------------------------------------------------
# CI column lookup: maps a variable name used on an axis to the column
# names holding its lower/upper 95% bound, where such a bound exists in
# the source data. Variables not listed here have no reported uncertainty
# (point estimates only) and are plotted without error bars, with a
# caption note added in the Rmd instead of fabricating one.
# -----------------------------------------------------------------------

ci_lookup <- list(
  "lambda_hat"                  = list(lwr = "lambda_lwr",               upr = "lambda_upr"),
  "Recombination probability"   = list(lwr = "asymp.LCL",                upr = "asymp.UCL"),
  "rate_genes1_median"          = list(lwr = "rate_genes1_lower_CI",     upr = "rate_genes1_upper_CI"),
  "prop_genes2_median"          = list(lwr = "prop_genes2_lower_CI",     upr = "prop_genes2_upper_CI"),
  "Y-intercept (j)_Estimate"    = list(lwr = "Y-intercept (j)_Lower CI", upr = "Y-intercept (j)_Upper 95% CI"),
  "yint_times_genome"           = list(lwr = "yint_times_genome_lwr",    upr = "yint_times_genome_upr")
)

# -----------------------------------------------------------------------
# fit_simple_lm(): a single, ungrouped linear model between two continuous
# variables, returning tidy statistics for annotation. This is always
# fitted on the FULL dataset passed in (never per-colour/per-shape
# subgroup), regardless of how points are subsequently coloured/shaped for
# display.
# -----------------------------------------------------------------------

fit_simple_lm <- function(df, x_var, y_var) {
  model_df <- df %>%
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))

  mod   <- stats::lm(stats::reformulate(x_var, y_var), data = model_df)
  glance_df <- broom::glance(mod)
  cor_test  <- stats::cor.test(model_df[[x_var]], model_df[[y_var]])

  list(
    model   = mod,
    n       = nrow(model_df),
    r       = unname(cor_test$estimate),
    r2      = glance_df$r.squared[1],
    p_value = glance_df$p.value[1]
  )
}

# -----------------------------------------------------------------------
# plot_correlation(): the main scatterplot workhorse.
#
# Key design choices (each fixes a specific earlier bug):
#   - colour/shape are mapped ONLY on the point/error-bar/label layers, not
#     in the top-level ggplot(aes(...)) call, and geom_smooth() is given an
#     explicit aes(group = 1). This guarantees a SINGLE, ungrouped
#     regression line/ribbon is drawn, regardless of how many
#     colour/shape categories are present - earlier versions implicitly
#     inherited colour/shape as grouping variables and silently fitted and
#     drew one regression per subgroup.
#   - the stats label is placed in dedicated head-room above all data
#     (expanding the upper y-limit and annotating at y = Inf), so it can
#     never overlap the fitted line or data points regardless of the sign
#     or strength of the correlation.
#   - geom_label_repel() is given an explicit colour aesthetic, since it no
#     longer inherits one from the top-level mapping.
#   - y-axis truncation (coord_cartesian) is optional and purely visual:
#     the regression is always fitted on the untruncated data first.
# -----------------------------------------------------------------------

plot_correlation <- function(df,
                             x_var, y_var,
                             x_lab = x_var, y_lab = y_var,
                             label_var   = "display_name",
                             colour_var  = "Transformability",
                             shape_var   = "Niche",
                             ylim        = NULL,
                             label_digits = 3) {

  plot_df <- df %>% dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))

  fit <- fit_simple_lm(plot_df, x_var, y_var)

  stats_lab <- safe_plotmath_label(list(
    plotmath_relation("italic(r)",      fit$r,       label_digits),
    plotmath_relation("italic(R)^2",    fit$r2,      label_digits),
    plotmath_relation("italic(p)",      fit$p_value, label_digits)
  ))

  has_x_ci <- x_var %in% names(ci_lookup) &&
    all(unlist(ci_lookup[[x_var]]) %in% names(plot_df))
  has_y_ci <- y_var %in% names(ci_lookup) &&
    all(unlist(ci_lookup[[y_var]]) %in% names(plot_df))

  p <- ggplot(plot_df, aes(x = .data[[x_var]], y = .data[[y_var]]))

  if (has_x_ci) {
    ci <- ci_lookup[[x_var]]
    p <- p + geom_errorbarh(
      aes(xmin = .data[[ci$lwr]], xmax = .data[[ci$upr]],
          colour = .data[[colour_var]]),
      height = 0, alpha = 0.5
    )
  }
  if (has_y_ci) {
    ci <- ci_lookup[[y_var]]
    p <- p + geom_errorbar(
      aes(ymin = .data[[ci$lwr]], ymax = .data[[ci$upr]],
          colour = .data[[colour_var]]),
      width = 0, alpha = 0.5
    )
  }

  p <- p +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
               colour = "black", fill = "grey60", linewidth = 0.8, alpha = 0.25) +
    geom_point(aes(colour = .data[[colour_var]], shape = .data[[shape_var]]),
              size = 2.5) +
    ggrepel::geom_label_repel(
      aes(label = .data[[label_var]], colour = .data[[colour_var]]),
      size = 3, parse = TRUE, force_pull = 0.5, segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = niche_palette) +
    scale_colour_manual(values = transformability_palette, aesthetics = c("colour", "fill")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(x = x_lab, y = y_lab) +
    theme_pangenome() +
    guides_pangenome()

  # When the y-axis is truncated for display (coord_cartesian only clips the
  # *view*; it does not change the underlying scale range), an annotation
  # placed at y = Inf would be clipped out of the visible window along with
  # the data it would otherwise sit above. In that case the label is placed
  # just inside the top of the truncated range instead.
  if (!is.null(ylim)) {
    ann_y     <- ylim[2]
    ann_vjust <- 1.3
    # clip stays "on" (the default): tall error bars/points that extend
    # beyond the truncated view (this is the whole point of truncating)
    # are cleanly cut off at the panel edge, rather than spilling out over
    # the legend/title as they would with clip = "off".
    p <- p + coord_cartesian(ylim = ylim)
  } else {
    ann_y     <- Inf
    ann_vjust <- 1.4
  }

  p <- p +
    annotate(
      "text", x = -Inf, y = ann_y,
      label = stats_lab$label, parse = stats_lab$parse,
      hjust = -0.05, vjust = ann_vjust, size = 3.3
    )

  list(plot = p, fit = fit)
}

# -----------------------------------------------------------------------
# plot_violin_comparison(): categorical (e.g. Transformability) vs
# continuous comparison with a Wilcoxon rank-sum test, styled to match
# nickjcroucher/recombination-conflict's plot_violin_kw().
# -----------------------------------------------------------------------

plot_violin_comparison <- function(df,
                                   x_var, y_var,
                                   x_lab = x_var, y_lab = y_var,
                                   colour_var = "Transformability",
                                   shape_var  = "Niche") {

  plot_df <- df %>% dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))

  wtest <- stats::wilcox.test(plot_df[[y_var]] ~ plot_df[[x_var]])

  p_lab <- safe_plotmath_label(list(
    plotmath_relation("italic(p)", wtest$p.value, 3)
  ), sep = "' '")
  p_lab$label <- paste0("'Wilcoxon test '*", p_lab$label)

  ggplot(plot_df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_violin(draw_quantiles = 0.5, alpha = 0.3, colour = "grey70") +
    geom_jitter(aes(colour = .data[[colour_var]], shape = .data[[shape_var]]),
               width = 0.15, height = 0, size = 2.2, alpha = 0.85) +
    scale_shape_manual(values = niche_palette) +
    scale_colour_manual(values = transformability_palette, aesthetics = c("colour", "fill")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = x_lab, y = y_lab) +
    annotate(
      "text", x = -Inf, y = Inf,
      label = p_lab$label, parse = p_lab$parse,
      hjust = -0.05, vjust = 1.4, size = 3.3
    ) +
    theme_pangenome() +
    guides_pangenome() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}
