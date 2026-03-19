################################################################################
################### Meta-analysis: Gut-related metabolites #####################
################################################################################

# Prepared by: Juliana Hernandez || j.a.hernandezvargas@umcutrecht.nl

################################################################################
### 1. Loading libraries and setting working directory
################################################################################

pacman::p_load(readxl, dplyr, stringr, purrr, metafor, ggplot2, janitor, writexl)

# Base directory (Paper 1_Metabolites)
base_dir <- "~/Documents/PhD files/Thesis/1. SR_GM and CVA/Analysis /Paper1_Metabolites"

# Excel file with harmonised data
xlsx_path <- file.path(base_dir, "Harmonisation for metaanalysis conversions.xlsx")
xlsx_pathrct <- file.path(base_dir, "Harmonisation RCTs.xlsx")

# Output directories
output_dir <- file.path(base_dir, "Meta_analysis_outputs")
plots_dir  <- file.path(output_dir, "Forest_plots")

# Create folders if they do not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

################################################################################
### 2. Load and filter metadata
################################################################################

#### Observational studies 
meta_xlsx <- readxl::read_excel(xlsx_path)

# Check available sheets
available_sheets <- readxl::excel_sheets(xlsx_path)
print(available_sheets)

# Sheet names to import
sheets_to_import <- c(
  "C-C_TMAO and CAD",
  "C-S_TMAO and CAD&burden",
  "C-S_TMAO and CAD extent plaque",
  "C-S_TMAO and CAD_all",
  "Cohorts_TMAO and cIMT",
  "Cohorts_TMAO and ASCVD",
  "C-S_TMAO and HTN",
  "C-C_TMAOpathway&CAD"
)

# Warn if any are missing (but still try to import the ones that exist)
missing_sheets <- setdiff(sheets_to_import, available_sheets)
if (length(missing_sheets) > 0) {
  warning("These sheets were not found in the Excel file:\n", paste(missing_sheets, collapse = "\n"))
}

# Import each sheet as a separate dataframe
cc_tmao_cad <- readxl::read_excel(xlsx_path, sheet = "C-C_TMAO and CAD")

cs_tmao_cad_all <- readxl::read_excel(xlsx_path, sheet = "C-S_TMAO and CAD_all")
                                  
cs_tmao_cad_burden <- readxl::read_excel(xlsx_path, sheet = "C-S_TMAO and CAD&burden")

cs_tmao_cad_extent <- readxl::read_excel(xlsx_path, sheet = "C-S_TMAO and CAD extent plaque")

cohorts_tmao_cimt <- readxl::read_excel(xlsx_path, sheet = "Cohorts_TMAO and cIMT")

cohorts_tmao_ascvd <- readxl::read_excel(xlsx_path, sheet = "Cohorts_TMAO and ASCVD")

cs_tmao_htn <- readxl::read_excel(xlsx_path, sheet = "C-S_TMAO and HTN")

cc_tmao_pathway_cad <- readxl::read_excel(xlsx_path, sheet = "C-C_TMAOpathway&CAD")


# Quick checks
dim(cc_tmao_cad)
dim (cs_tmao_cad_all)
dim(cs_tmao_cad_burden)
dim(cs_tmao_cad_extent)
dim(cohorts_tmao_cimt)
dim(cohorts_tmao_ascvd)
dim(cs_tmao_htn)
dim(cc_tmao_pathway_cad)

#### RCTs
meta_rcts <- readxl::read_excel(xlsx_pathrct)

################################################################################
### 3. Running meta-analyses
################################################################################

################################################################################
# 3.1. Forest plot for TMAO and CAD (Case-Control studies) — PER 1 SD (TMAO)
################################################################################

#--------------------------------
# Helper functions
#--------------------------------

# CI formatter
fmt_ci <- function(or, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", or, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

# ---- Load data ----
df_raw <- read_excel(xlsx_path, sheet = "C-C_TMAO and CAD") %>%
  clean_names()

# ---- Prepare dataset (per-1SD columns) ----
df <- df_raw %>%
  mutate(
    OR       = as.numeric(or_per1sd),
    LCL      = as.numeric(ci_low_per1sd),
    UCL      = as.numeric(ci_high_per1sd),
    N        = suppressWarnings(as.integer(sample_size)),
    Cases    = suppressWarnings(as.integer(number_cases)),
    Controls = suppressWarnings(as.integer(number_controls)),
    Study    = paste0(first_author, " (", year, ")"),
    yi       = log(OR),
    sei      = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi       = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

# ---- Random-effects meta-analysis ----
res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- Meta package object for reporting ----
m_re <- meta::metagen(
  TE      = df$yi,                 # log(OR)
  seTE    = sqrt(df$vi),
  studlab = df$Study,
  sm      = "OR",
  comb.fixed  = FALSE,
  comb.random = TRUE,
  method.tau  = "REML",
  hakn        = FALSE
)

# ---- Heterogeneity (from metafor) ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100
tau2 <- as.numeric(res_re$tau2)

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "; I²=", sprintf("%.1f", I2), "%; Tau²=", sprintf("%.3f", tau2)
)

# ---- Random-effects weights ----
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

################################################################################
# Forest plot layout 
################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- x-axis ----
x_ticks <- c(0.75, 1, 1.25, 1.5, 2, 3, 4)
x_at    <- log(x_ticks)
alim    <- log(c(0.75, 4.0))

# ---- Expand forest region ----
xlim <- c(alim[1] - 2.6, alim[2] + 2.1)   # was -3.8 / +3.3
ylim <- c(-1.7, k + 2.2)

# ---- Column positions ----
x_study    <- xlim[1] + 0.15
x_n        <- alim[1] - 0.95
x_casectrl <- alim[1] - 0.20
x_weight   <- alim[2] + 0.35
x_est      <- alim[2] + 1.03

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

# ---- Header line position ----
y_header_line <- k + 1.1

################################################################################
# Plot function
################################################################################

make_forestplot <- function() {
  
  par(mar = c(3.8, 6, 2.5, 8))
  par(xpd = NA)
  
  # Base forest without refline
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at,
    atransf = exp,
    refline = NA,
    xlab    = "",
    rows    = rows,
    psize   = 1.2,
    cex     = 1.1,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE
  )
  
  usr <- par("usr")
  
  # ---- Alternating row shading ----
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(usr[1], rows[i] - 0.5, usr[2], rows[i] + 0.5,
           col = bg_alt, border = NA)
    }
  }
  
  # ---- Clean header band ----
  rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # ---- Null line at OR=1  ----
  segments(0, usr[3], 0, y_header_line, lty = 2)
  
  # ---- CI lines ----
  segments(df$yi - 1.96 * sqrt(df$vi), rows,
           df$yi + 1.96 * sqrt(df$vi), rows,
           col = col_ci, lwd = 2.2)
  
  segments(df$yi - 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi - 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  segments(df$yi + 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi + 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  # ---- Squares sized by RE weights ----
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$yi, rows, pch = 15, cex = 1.0 + 1.4 * w_scaled, col = col_sq)
  
  # ---- Column headers ----
  text(x_study,    k + 1.5, "Study",                font = 2, adj = 0)
  text(x_n,        k + 1.5, "N",                    font = 2, adj = 1)
  text(x_casectrl, k + 1.5, "Cases/Controls",       font = 2, adj = 1)
  text(x_weight,   k + 1.5, "Weight (%)",           font = 2, adj = 0)
  text(x_est,      k + 1.5, "OR per 1 SD (95% CI)", font = 2, adj = 0)
  
  # ---- Row labels ----
  text(x_study, rows, df$Study, adj = 0)
  text(x_n, rows, df$N, adj = 1)
  text(x_casectrl, rows, paste0(df$Cases, " / ", df$Controls), adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est, rows, fmt_ci(df$OR, df$LCL, df$UCL), adj = 0)
  
  # ---- Random-effects diamond ----
  metafor::addpoly(
    res_re,
    row = row_re,
    atransf = exp,
    mlab = "",
    col = col_dia,
    border = col_dia
  )
  
  # ---- Random-effects label and estimate ----
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- exp(c(res_re$b, res_re$ci.lb, res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "—", adj = 0, cex = 1.10)
  
  # ---- Heterogeneity ----
  text(x_study, row_het, het_txt, adj = 0, cex = 1.07)
}

################################################################################
# Save outputs
################################################################################

pdf(file.path(plots_dir, "Forestplot_TMAO_CAD_per1SD.pdf"), width = 14, height = 9)
make_forestplot()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_TMAO_CAD_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_TMAO_CAD_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot()
dev.off()

message("Saved forest plots (PDF, TIFF, JPEG) to: ", plots_dir)

# Summary report 
print(summary(m_re))

################################################################################
# Cochrane Diagnostics: Sensitivity & Publication Bias
################################################################################

################################################################################
### Publication-Ready Diagnostics (Sensitivity & Bias)
################################################################################

# Create a subdirectory for diagnostics 
diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo <- metafor::leave1out(res_re)

if (!"pval" %in% names(res_loo)) {
  res_loo$pval <- 2 * pnorm(abs(res_loo$zval), lower.tail = FALSE)
}

loo_est  <- exp(res_loo$estimate)
loo_lb   <- exp(res_loo$ci.lb)
loo_ub   <- exp(res_loo$ci.ub)
loo_pval <- res_loo$pval
studies  <- df$Study

k <- res_re$k
y_pos <- rev(seq_len(k))

forest_xmin  <- min(loo_lb) * 0.92
forest_xmax  <- max(loo_ub) * 1.08
forest_range <- forest_xmax - forest_xmin

study_col_x  <- forest_xmin - 1.05 * forest_range
forest_left  <- forest_xmin
forest_right <- forest_xmax
effect_col_x <- forest_right + 0.15 * forest_range
p_col_x      <- forest_right + 0.82 * forest_range
plot_xmax    <- forest_right + 1.20 * forest_range

png(file.path(diag_dir, "Sensitivity_LOO_TMAO_CAD.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 4, 1.2, 1.2), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled OR if study is excluded",
     ylab = "",
     bty = "n")

axis_ticks <- pretty(c(forest_left, forest_right))
axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)

text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

segments(loo_lb, y_pos, loo_ub, y_pos, col = "#19C3A3", lwd = 1.4)
points(loo_est, y_pos, pch = 16, cex = 0.9, col = "#19C3A3")

segments(x0 = exp(res_re$b), y0 = 0.5, x1 = exp(res_re$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

text(x = study_col_x, y = k + 0.72, labels = "Omitted study", adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "OR with 95% CI", adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value", adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot  ----
yi  <- as.numeric(res_re$yi)
sei <- sqrt(res_re$vi)
mu  <- as.numeric(res_re$b)

# Change this label depending on the analysis:
# xlab_funnel <- "Mean difference"
xlab_funnel <- "Log Odds Ratio"

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq), color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = xlab_funnel,
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_TMAO_CAD.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
cat("\n--- Publication Bias Statistical Tests ---\n")

egger_result <- metafor::regtest(res_re, model = "rma", predictor = "sei")
print(egger_result)

begg_result <- metafor::ranktest(res_re)
print(begg_result)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"))
cat("Egger's Regression Test:\n")
print(egger_result)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_result)
sink()

message("Sensitivity plot, Funnel plot, and statistical tests saved to: ", diag_dir)

#############################################################################################
# 3.2a. Forest plot for TMAO PATHWAY (betaine + phenylacetylglutamine) and CAD (Case-control)
#############################################################################################

#--------------------------------
# Helper functions
#--------------------------------
fmt_ci <- function(or, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", or, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

df_raw <- read_excel(xlsx_path, sheet = "C-C_TMAOpathway&CAD") %>%
  clean_names() %>%
  dplyr::filter(
    stringr::str_to_lower(gut_metabolite) %in% c(
      "betaine",
      "phenylacetylglutamine",
      "phenylacethylglutamine"
    )
  )

df <- df_raw %>%
  mutate(
    OR       = as.numeric(or_per1sd),
    LCL      = as.numeric(ci_low_per1sd),
    UCL      = as.numeric(ci_high_per1sd),
    N        = suppressWarnings(as.integer(sample_size)),
    Cases    = suppressWarnings(as.integer(number_cases)),
    Controls = suppressWarnings(as.integer(number_controls)),
    Study    = paste0(first_author, " (", year, ")"),
    yi       = log(OR),
    sei      = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi       = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

# ---- Random-effects meta-analysis  ----
res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- meta package object for reporting ----
m_re <- meta::metagen(
  TE      = df$yi,                 # log(OR)
  seTE    = sqrt(df$vi),
  studlab = df$Study,
  sm      = "OR",
  comb.fixed  = FALSE,
  comb.random = TRUE,
  method.tau  = "REML",
  hakn        = FALSE
)

# ---- Heterogeneity (from metafor) ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100
tau2 <- as.numeric(res_re$tau2)

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "; I²=", sprintf("%.1f", I2), "%; Tau²=", sprintf("%.3f", tau2)
)

# ---- Random-effects weights ----
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

################################################################################
# Forest plot layout
################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- x-axis----
x_ticks <- c(0.50, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)
x_at    <- log(x_ticks)
alim    <- log(c(0.50, 3.0))

# ---- Fixing forest region ----
xlim <- c(alim[1] - 3.2, alim[2] + 2.1)
ylim <- c(-1.7, k + 2.2)

# ---- Column positions ----
x_study    <- xlim[1] + 0.15
x_n        <- alim[1] - 1.25   
x_casectrl <- alim[1] - 0.20
x_weight   <- alim[2] + 0.35
x_est      <- alim[2] + 1.03

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

# ---- Header line position ----
y_header_line <- k + 1.1

################################################################################
# Plot function
################################################################################

make_forestplot_pathway <- function() {
  
  par(mar = c(3.8, 6, 2.5, 8))
  par(xpd = NA)
  
  # Base forest without refline
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at,
    atransf = exp,
    refline = NA,
    xlab    = "",
    rows    = rows,
    psize   = 1.2,
    cex     = 1.1,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE
  )
  
  usr <- par("usr")
  
  # ---- Alternating row shading ----
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(usr[1], rows[i] - 0.5, usr[2], rows[i] + 0.5,
           col = bg_alt, border = NA)
    }
  }
  
  # ---- Clean header band ----
  rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # ---- Null line at OR=1 ----
  segments(0, usr[3], 0, y_header_line, lty = 2)
  
  # ---- CI lines + caps ----
  segments(df$yi - 1.96 * sqrt(df$vi), rows,
           df$yi + 1.96 * sqrt(df$vi), rows,
           col = col_ci, lwd = 2.2)
  
  segments(df$yi - 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi - 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  segments(df$yi + 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi + 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  # ---- Squares sized by RE weights ----
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$yi, rows, pch = 15, cex = 1.0 + 1.4 * w_scaled, col = col_sq)
  
  # ---- Column headers ----
  text(x_study,    k + 1.5, "Study",                font = 2, adj = 0)
  text(x_n,        k + 1.5, "N",                    font = 2, adj = 1)
  text(x_casectrl, k + 1.5, "Cases/Controls",       font = 2, adj = 1)
  text(x_weight,   k + 1.5, "Weight (%)",           font = 2, adj = 0)
  text(x_est,      k + 1.5, "OR per 1 SD (95% CI)", font = 2, adj = 0)
  
  # ---- Row labels ----
  text(x_study, rows, df$Study, adj = 0)
  text(x_n, rows, df$N, adj = 1)
  text(x_casectrl, rows, paste0(df$Cases, " / ", df$Controls), adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est, rows, fmt_ci(df$OR, df$LCL, df$UCL), adj = 0)
  
  # ---- Random-effects diamond ----
  metafor::addpoly(
    res_re,
    row = row_re,
    atransf = exp,
    mlab = "",
    col = col_dia,
    border = col_dia
  )
  
  # ---- Random-effects label + estimate ----
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- exp(c(res_re$b, res_re$ci.lb, res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "—", adj = 0, cex = 1.10)
  
  # ---- Heterogeneity ----
  text(x_study, row_het, het_txt, adj = 0, cex = 1.07)
}

################################################################################
# Save outputs
################################################################################

pdf(file.path(plots_dir, "Forestplot_TMAO_pathway_CAD_per1SD.pdf"), width = 14, height = 9)
make_forestplot_pathway()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_TMAO_pathway_CAD_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot_pathway()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_TMAO_pathway_CAD_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot_pathway()
dev.off()

message("Saved forest plots (PDF, TIFF, JPEG) to: ", plots_dir)

# --- Printing summary ---
print(summary(m_re))

################################################################################
### Publication-Ready Diagnostics: TMAO Pathway & CAD
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- Rebuild the pathway-specific dataset and model  ----
df_pathway <- read_excel(xlsx_path, sheet = "C-C_TMAOpathway&CAD") %>%
  clean_names() %>%
  dplyr::filter(
    stringr::str_to_lower(gut_metabolite) %in% c(
      "betaine",
      "phenylacetylglutamine",
      "phenylacethylglutamine"
    )
  ) %>%
  mutate(
    OR       = as.numeric(or_per1sd),
    LCL      = as.numeric(ci_low_per1sd),
    UCL      = as.numeric(ci_high_per1sd),
    N        = suppressWarnings(as.integer(sample_size)),
    Cases    = suppressWarnings(as.integer(number_cases)),
    Controls = suppressWarnings(as.integer(number_controls)),
    Study    = paste0(first_author, " (", year, ")"),
    yi       = log(OR),
    sei      = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi       = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df_pathway) >= 2)

res_re_pathway <- metafor::rma(yi, vi, data = df_pathway, method = "REML")

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo_pathway <- metafor::leave1out(res_re_pathway)

if (!"pval" %in% names(res_loo_pathway)) {
  res_loo_pathway$pval <- 2 * pnorm(abs(res_loo_pathway$zval), lower.tail = FALSE)
}

loo_est  <- exp(res_loo_pathway$estimate)
loo_lb   <- exp(res_loo_pathway$ci.lb)
loo_ub   <- exp(res_loo_pathway$ci.ub)
loo_pval <- res_loo_pathway$pval
studies  <- df_pathway$Study

k <- nrow(df_pathway)
y_pos <- rev(seq_len(k))

forest_xmin  <- min(loo_lb) * 0.92
forest_xmax  <- max(loo_ub) * 1.08
forest_range <- forest_xmax - forest_xmin

study_col_x  <- forest_xmin - 1.05 * forest_range
forest_left  <- forest_xmin
forest_right <- forest_xmax
effect_col_x <- forest_right + 0.15 * forest_range
p_col_x      <- forest_right + 0.82 * forest_range
plot_xmax    <- forest_right + 1.20 * forest_range

png(file.path(diag_dir, "Sensitivity_LOO_TMAO_Pathway_CAD.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 4, 1.2, 1.2), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled OR if study is excluded",
     ylab = "",
     bty = "n")

axis_ticks <- pretty(c(forest_left, forest_right))
axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)

text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

segments(loo_lb, y_pos, loo_ub, y_pos, col = "#19C3A3", lwd = 1.4)
points(loo_est, y_pos, pch = 16, cex = 0.9, col = "#19C3A3")

segments(x0 = exp(res_re_pathway$b), y0 = 0.5, x1 = exp(res_re_pathway$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
     adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "OR with 95% CI",
     adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value",
     adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot ----
yi  <- as.numeric(res_re_pathway$yi)
sei <- sqrt(res_re_pathway$vi)
mu  <- as.numeric(res_re_pathway$b)

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
        color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = "Log Odds Ratio",
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_TMAO_Pathway_CAD.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
egger_pathway <- metafor::regtest(res_re_pathway, model = "rma", predictor = "sei")
begg_pathway  <- metafor::ranktest(res_re_pathway)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
cat("\n\n********************************************************\n")
cat("ANALYSIS: TMAO Pathway (Betaine + PAG) and CAD\n")
cat("********************************************************\n")
cat("Number of studies included: ", nrow(df_pathway), "\n", sep = "")
cat("Studies included: ", paste(df_pathway$Study, collapse = "; "), "\n\n", sep = "")
cat("Egger's Regression Test:\n")
print(egger_pathway)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_pathway)
sink()

message("TMAO Pathway diagnostics saved to the 'Diagnostics' folder.")

#############################################################################################
# 3.2b. Forest plot for TMAO PATHWAY (choline + phenylacetylglutamine) and CAD (Case-control)
#############################################################################################

#--------------------------------
# Helper functions
#--------------------------------
fmt_ci <- function(or, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", or, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

# ---- Load data ----
df_raw <- read_excel(xlsx_path, sheet = "C-C_TMAOpathway&CAD") %>%
  clean_names() %>%
  dplyr::filter(
    stringr::str_to_lower(gut_metabolite) %in% c(
      "choline",
      "phenylacetylglutamine",
      "phenylacethylglutamine"
    )
  )

# ---- Prepare dataset (per-1SD columns) ----
df <- df_raw %>%
  mutate(
    OR       = as.numeric(or_per1sd),
    LCL      = as.numeric(ci_low_per1sd),
    UCL      = as.numeric(ci_high_per1sd),
    N        = suppressWarnings(as.integer(sample_size)),
    Cases    = suppressWarnings(as.integer(number_cases)),
    Controls = suppressWarnings(as.integer(number_controls)),
    Study    = paste0(first_author, " (", year, ")"),
    yi       = log(OR),
    sei      = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi       = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

# ---- Random-effects meta-analysis (metafor) ----
res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- Meta package object for reporting ----
m_re <- meta::metagen(
  TE      = df$yi,                 # log(OR)
  seTE    = sqrt(df$vi),
  studlab = df$Study,
  sm      = "OR",
  comb.fixed  = FALSE,
  comb.random = TRUE,
  method.tau  = "REML",
  hakn        = FALSE
)

# ---- Heterogeneity ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100
tau2 <- as.numeric(res_re$tau2)

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "; I²=", sprintf("%.1f", I2), "%; Tau²=", sprintf("%.3f", tau2)
)

# ---- Random-effects weights ----
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

################################################################################
# Forest plot layout 
################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- Fixed x-axis  ----
x_ticks <- c(0.30, 0.40, 0.50, 0.60, 0.75, 1, 1.25, 1.5, 2, 2.5)
x_at    <- log(x_ticks)
alim    <- log(c(0.30, 2.5))

# ---- Keep the "shift forest right" behavior ----
xlim <- c(alim[1] - 3.2, alim[2] + 2.1)
ylim <- c(-1.7, k + 2.2)

# ---- Column positions ----
x_study    <- xlim[1] + 0.15
x_n        <- alim[1] - 1.25
x_casectrl <- alim[1] - 0.20
x_weight   <- alim[2] + 0.35
x_est      <- alim[2] + 1.03

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

# ---- Header line position ----
y_header_line <- k + 1.1

################################################################################
# Plot function
################################################################################

make_forestplot_pathway2 <- function() {
  
  par(mar = c(3.8, 6, 2.5, 8))
  par(xpd = NA)
  
  # Base forest without refline
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at,
    atransf = exp,
    refline = NA,
    xlab    = "",
    rows    = rows,
    psize   = 1.2,
    cex     = 1.1,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE
  )
  
  usr <- par("usr")
  
  # ---- Zebra shading ----
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(usr[1], rows[i] - 0.5, usr[2], rows[i] + 0.5,
           col = bg_alt, border = NA)
    }
  }
  
  # ---- Clean header band ----
  rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # ---- Null line at OR=1  ----
  segments(0, usr[3], 0, y_header_line, lty = 2)
  
  # ---- CI lines + caps ----
  segments(df$yi - 1.96 * sqrt(df$vi), rows,
           df$yi + 1.96 * sqrt(df$vi), rows,
           col = col_ci, lwd = 2.2)
  segments(df$yi - 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi - 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  segments(df$yi + 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi + 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  # ---- Squares ----
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$yi, rows, pch = 15, cex = 1.0 + 1.4 * w_scaled, col = col_sq)
  
  # ---- Column headers ----
  text(x_study,    k + 1.5, "Study",                font = 2, adj = 0)
  text(x_n,        k + 1.5, "N",                    font = 2, adj = 1)
  text(x_casectrl, k + 1.5, "Cases/Controls",       font = 2, adj = 1)
  text(x_weight,   k + 1.5, "Weight (%)",           font = 2, adj = 0)
  text(x_est,      k + 1.5, "OR per 1 SD (95% CI)", font = 2, adj = 0)
  
  # ---- Row labels ----
  text(x_study, rows, df$Study, adj = 0)
  text(x_n, rows, df$N, adj = 1)
  text(x_casectrl, rows, paste0(df$Cases, " / ", df$Controls), adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est, rows, fmt_ci(df$OR, df$LCL, df$UCL), adj = 0)
  
  # ---- Random-effects diamond ----
  metafor::addpoly(
    res_re,
    row = row_re,
    atransf = exp,
    mlab = "",
    col = col_dia,
    border = col_dia
  )
  
  # ---- Random-effects label + estimate ----
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- exp(c(res_re$b, res_re$ci.lb, res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "—", adj = 0, cex = 1.10)
  
  # ---- Heterogeneity ----
  text(x_study, row_het, het_txt, adj = 0, cex = 1.07)
}

##################################################################################
# Save in multiple formats
##################################################################################

pdf(file.path(plots_dir, "Forestplot_TMAO_pathway_2_CAD_per1SD.pdf"), width = 14, height = 9)
make_forestplot_pathway2()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_TMAO_pathway_2_CAD_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot_pathway2()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_TMAO_pathway_2_CAD_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot_pathway2()
dev.off()

message("Saved TMAO pathway 2 forest plots (PDF, TIFF, JPEG) to: ", plots_dir)

# --- Printing summary ---
print(summary(m_re))

################################################################################
### Publication-Ready Diagnostics: TMAO Pathway 2 (Choline + PAG) & CAD
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- Rebuild the pathway-2-specific dataset and model  ----
df_pathway2 <- read_excel(xlsx_path, sheet = "C-C_TMAOpathway&CAD") %>%
  clean_names() %>%
  dplyr::filter(
    stringr::str_to_lower(gut_metabolite) %in% c(
      "choline",
      "phenylacetylglutamine",
      "phenylacethylglutamine"
    )
  ) %>%
  mutate(
    OR       = as.numeric(or_per1sd),
    LCL      = as.numeric(ci_low_per1sd),
    UCL      = as.numeric(ci_high_per1sd),
    N        = suppressWarnings(as.integer(sample_size)),
    Cases    = suppressWarnings(as.integer(number_cases)),
    Controls = suppressWarnings(as.integer(number_controls)),
    Study    = paste0(first_author, " (", year, ")"),
    yi       = log(OR),
    sei      = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi       = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df_pathway2) >= 2)

res_re_pathway2 <- metafor::rma(yi, vi, data = df_pathway2, method = "REML")

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo_pathway2 <- metafor::leave1out(res_re_pathway2)

if (!"pval" %in% names(res_loo_pathway2)) {
  res_loo_pathway2$pval <- 2 * pnorm(abs(res_loo_pathway2$zval), lower.tail = FALSE)
}

loo_est  <- exp(res_loo_pathway2$estimate)
loo_lb   <- exp(res_loo_pathway2$ci.lb)
loo_ub   <- exp(res_loo_pathway2$ci.ub)
loo_pval <- res_loo_pathway2$pval
studies  <- df_pathway2$Study

k <- nrow(df_pathway2)
y_pos <- rev(seq_len(k))

forest_xmin  <- min(loo_lb) * 0.92
forest_xmax  <- max(loo_ub) * 1.08
forest_range <- forest_xmax - forest_xmin

study_col_x  <- forest_xmin - 1.05 * forest_range
forest_left  <- forest_xmin
forest_right <- forest_xmax
effect_col_x <- forest_right + 0.15 * forest_range
p_col_x      <- forest_right + 0.82 * forest_range
plot_xmax    <- forest_right + 1.20 * forest_range

png(file.path(diag_dir, "Sensitivity_LOO_TMAO_Pathway2_CAD.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 4, 1.2, 1.2), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled OR if study is excluded",
     ylab = "",
     bty = "n")

axis_ticks <- pretty(c(forest_left, forest_right))
axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)

text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

segments(loo_lb, y_pos, loo_ub, y_pos, col = "#19C3A3", lwd = 1.4)
points(loo_est, y_pos, pch = 16, cex = 0.9, col = "#19C3A3")

segments(x0 = exp(res_re_pathway2$b), y0 = 0.5, x1 = exp(res_re_pathway2$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
     adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "OR with 95% CI",
     adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value",
     adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot ----
yi  <- as.numeric(res_re_pathway2$yi)
sei <- sqrt(res_re_pathway2$vi)
mu  <- as.numeric(res_re_pathway2$b)

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
        color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = "Log Odds Ratio",
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_TMAO_Pathway2_CAD.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
egger_pathway2 <- metafor::regtest(res_re_pathway2, model = "rma", predictor = "sei")
begg_pathway2  <- metafor::ranktest(res_re_pathway2)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
cat("\n\n********************************************************\n")
cat("ANALYSIS: TMAO Pathway 2 (Choline + PAG) and CAD\n")
cat("********************************************************\n")
cat("Number of studies included: ", nrow(df_pathway2), "\n", sep = "")
cat("Studies included: ", paste(df_pathway2$Study, collapse = "; "), "\n\n", sep = "")
cat("Egger's Regression Test:\n")
print(egger_pathway2)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_pathway2)
sink()

message("Diagnostics for TMAO Pathway 2 saved to: ", diag_dir)

##################################################################################
# 3.3. Forest plot for COHORTS — TMAO and ASCVD (RR)
##################################################################################

#--------------------------------
# Helper functions
#--------------------------------
fmt_ci <- function(rr, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", rr, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

# ---- Load data ----
df_raw <- read_excel(xlsx_path, sheet = "Cohorts_TMAO and ASCVD") %>%
  clean_names()

# ---- Prepare dataset  ----
df <- df_raw %>%
  mutate(
    RR    = as.numeric(rr_per1sd),
    LCL   = as.numeric(ci_low_per1sd),
    UCL   = as.numeric(ci_high_per1sd),
    N     = suppressWarnings(as.integer(sample_size)),
    Study    = paste0(first_author, " (", year, ")"),
    yi    = log(RR),
    sei   = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi    = sei^2
  ) %>%
  filter(
    is.finite(RR), is.finite(LCL), is.finite(UCL),
    RR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

# ---- Random-effects meta-analysis (metafor) ----
res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- Meta package object for reporting ----
m_re <- meta::metagen(
  TE      = df$yi,
  seTE    = sqrt(df$vi),
  studlab = df$Study,
  sm      = "RR",
  comb.fixed  = FALSE,
  comb.random = TRUE,
  method.tau  = "REML",
  hakn        = FALSE
)

# ---- Heterogeneity ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100
tau2 <- as.numeric(res_re$tau2)

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "; I²=", sprintf("%.1f", I2), "%; Tau²=", sprintf("%.3f", tau2)
)

# ---- Random-effects weights ----
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

##################################################################################
# Layout
##################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- Axis ----
x_ticks_major <- seq(0.8, 2.4, by = 0.4)
x_ticks_minor <- seq(0.8, 2.4, by = 0.2)

x_at_major <- log(x_ticks_major)
x_at_minor <- log(x_ticks_minor)

alim <- log(c(0.75, 2.5))

# ---- Forest region ----
xlim <- c(alim[1] - 4.8, alim[2] + 3.8)
ylim <- c(-1.7, k + 2.2)

# ---- Column positions ----
x_study  <- xlim[1] + 0.50
x_n      <- alim[1] - 0.80
x_weight <- alim[2] + 0.85
x_est    <- alim[2] + 2.00

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

# ---- Header line position ----
y_header_line <- k + 1.1

##################################################################################
# Plot function
##################################################################################

make_forestplot_cohorts <- function() {
  
  par(mar = c(3.8, 6, 2.5, 9))
  par(xpd = NA)
  
  # Base forest without refline
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at_major,
    atransf = exp,
    refline = NA,
    xlab    = "",
    rows    = rows,
    psize   = 1.2,
    cex     = 1.1,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE
  )
  
  usr <- par("usr")
  
  # ---- Zebra shading ----
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(usr[1], rows[i] - 0.5, usr[2], rows[i] + 0.5,
           col = bg_alt, border = NA)
    }
  }
  
  # ---- Clean header band ----
  rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # ---- Null line at RR=1 ----
  segments(0, usr[3], 0, y_header_line, lty = 2)
  
  # ---- Minor ticks  ----
  axis(
    side   = 1,
    at     = x_at_minor,
    labels = FALSE,
    tcl    = -0.3
  )
  
  # ---- CI lines + caps ----
  segments(df$yi - 1.96 * sqrt(df$vi), rows,
           df$yi + 1.96 * sqrt(df$vi), rows,
           col = col_ci, lwd = 2.2)
  segments(df$yi - 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi - 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  segments(df$yi + 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi + 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  # ---- Squares sized by RE weights ----
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$yi, rows, pch = 15, cex = 1.0 + 1.4 * w_scaled, col = col_sq)
  
  # ---- Column headers ----
  text(x_study,  k + 1.5, "Study",                font = 2, adj = 0)
  text(x_n,      k + 1.5, "N",                    font = 2, adj = 1)
  text(x_weight, k + 1.5, "Weight (%)",           font = 2, adj = 0)
  text(x_est,    k + 1.5, "RR per 1 SD (95% CI)", font = 2, adj = 0)
  
  # ---- Row labels ----
  text(x_study, rows, df$Study, adj = 0)
  text(x_n, rows, df$N, adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est, rows, fmt_ci(df$RR, df$LCL, df$UCL), adj = 0)
  
  # ---- Random-effects diamond ----
  metafor::addpoly(
    res_re,
    row = row_re,
    atransf = exp,
    mlab = "",
    col = col_dia,
    border = col_dia
  )
  
  # ---- Random-effects label + estimate ----
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- exp(c(res_re$b, res_re$ci.lb, res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "—", adj = 0, cex = 1.10)
  
  # ---- Heterogeneity ----
  text(x_study, row_het, het_txt, adj = 0, cex = 1.07)
}

##################################################################################
# Save in multiple formats
##################################################################################

pdf(file.path(plots_dir, "Forestplot_COHORTS_TMAO_ASCVD_per1SD.pdf"), width = 14, height = 9)
make_forestplot_cohorts()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_COHORTS_TMAO_ASCVD_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot_cohorts()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_COHORTS_TMAO_ASCVD_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot_cohorts()
dev.off()

message("Saved COHORT TMAO–ASCVD forest plots (PDF, TIFF, JPEG) to: ", plots_dir)

# --- Printing summary ---
print(summary(m_re))

################################################################################
### Publication-Ready Diagnostics: Cohorts (TMAO and ASCVD)
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- Rebuild the cohort-specific dataset and model  ----
df_cohorts <- read_excel(xlsx_path, sheet = "Cohorts_TMAO and ASCVD") %>%
  clean_names() %>%
  mutate(
    RR    = as.numeric(rr_per1sd),
    LCL   = as.numeric(ci_low_per1sd),
    UCL   = as.numeric(ci_high_per1sd),
    N     = suppressWarnings(as.integer(sample_size)),
    Study = paste0(first_author, " (", year, ")"),
    yi    = log(RR),
    sei   = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi    = sei^2
  ) %>%
  filter(
    is.finite(RR), is.finite(LCL), is.finite(UCL),
    RR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df_cohorts) >= 2)

res_re_cohorts <- metafor::rma(yi, vi, data = df_cohorts, method = "REML")

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo_cohorts <- metafor::leave1out(res_re_cohorts)

if (!"pval" %in% names(res_loo_cohorts)) {
  res_loo_cohorts$pval <- 2 * pnorm(abs(res_loo_cohorts$zval), lower.tail = FALSE)
}

loo_est  <- exp(res_loo_cohorts$estimate)
loo_lb   <- exp(res_loo_cohorts$ci.lb)
loo_ub   <- exp(res_loo_cohorts$ci.ub)
loo_pval <- res_loo_cohorts$pval
studies  <- df_cohorts$Study

k <- nrow(df_cohorts)
y_pos <- rev(seq_len(k))

forest_xmin  <- min(loo_lb) * 0.92
forest_xmax  <- max(loo_ub) * 1.08
forest_range <- forest_xmax - forest_xmin

study_col_x  <- forest_xmin - 1.05 * forest_range
forest_left  <- forest_xmin
forest_right <- forest_xmax
effect_col_x <- forest_right + 0.15 * forest_range
p_col_x      <- forest_right + 0.82 * forest_range
plot_xmax    <- forest_right + 1.20 * forest_range

png(file.path(diag_dir, "Sensitivity_LOO_Cohorts_TMAO_ASCVD.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 4, 1.2, 1.2), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled RR if study is excluded",
     ylab = "",
     bty = "n")

axis_ticks <- pretty(c(forest_left, forest_right))
axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)

text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

segments(loo_lb, y_pos, loo_ub, y_pos, col = "#19C3A3", lwd = 1.4)
points(loo_est, y_pos, pch = 16, cex = 0.9, col = "#19C3A3")

segments(x0 = exp(res_re_cohorts$b), y0 = 0.5, x1 = exp(res_re_cohorts$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
     adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "RR with 95% CI",
     adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value",
     adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot ----
yi  <- as.numeric(res_re_cohorts$yi)
sei <- sqrt(res_re_cohorts$vi)
mu  <- as.numeric(res_re_cohorts$b)

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
        color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = "Log Risk Ratio",
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_Cohorts_TMAO_ASCVD.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
egger_cohorts <- metafor::regtest(res_re_cohorts, model = "rma", predictor = "sei")
begg_cohorts  <- metafor::ranktest(res_re_cohorts)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
cat("\n\n********************************************************\n")
cat("ANALYSIS: COHORTS - TMAO and ASCVD (RR)\n")
cat("********************************************************\n")
cat("Number of studies included: ", nrow(df_cohorts), "\n", sep = "")
cat("Studies included: ", paste(df_cohorts$Study, collapse = "; "), "\n\n", sep = "")
cat("Egger's Regression Test:\n")
print(egger_cohorts)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_cohorts)
sink()

message("Diagnostics for Cohorts saved to: ", diag_dir)

##################################################################################
# 3.4. Forest plot for CROSS-SECTIONAL — TMAO and HTN (OR)
##################################################################################

#--------------------------------
# Helper functions
#--------------------------------
fmt_ci <- function(or, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", or, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

# ---- Load data ----
df_raw <- read_excel(xlsx_path, sheet = "C-S_TMAO and HTN") %>%
  clean_names()

# ---- Prepare dataset ----
df <- df_raw %>%
  mutate(
    OR  = suppressWarnings(as.numeric(dplyr::coalesce(or_per1sd, or))),
    LCL = suppressWarnings(as.numeric(dplyr::coalesce(ci_low_per1sd, or_ci_lower))),
    UCL = suppressWarnings(as.numeric(dplyr::coalesce(ci_high_per1sd, or_ci_upper))),
    N   = suppressWarnings(as.integer(sample_size)),
    Study    = paste0(first_author, " (", year, ")"),
    yi  = log(OR),
    sei = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi  = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

# ---- Random-effects meta-analysis ----
res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- Weights (random-effects) ----
tau2 <- as.numeric(res_re$tau2)
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

# ---- Heterogeneity ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "\nI\u00B2=", sprintf("%.1f", I2), "%; Tau\u00B2=", sprintf("%.3f", tau2)
)

##################################################################################
# Layout
##################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- OR grid: 0.2 to 100 ----
x_ticks <- c(0.2, 0.25, 0.5, 1, 2, 5, 10, 25, 50, 100)
x_at    <- log(x_ticks)
alim    <- log(c(0.2, 100))

# ---- Plot region ----
ylim <- c(-1.3, k + 1.1)
xlim <- c(alim[1] - 10.5, alim[2] + 10.0)

# ---- Column positions ----
x_study  <- xlim[1] + 0.50
x_n      <- alim[1] - 0.60  
x_weight <- alim[2] + 3.00
x_est    <- alim[2] + 5.40

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

# ---- Header line position ----
y_header_line <- max(rows) + 0.25

##################################################################################
# Plot function
##################################################################################

make_forestplot_crosssectional <- function() {
  
  par(mar = c(3.8, 6, 2.5, 8))
  par(xpd = NA)
  
  # Base forest 
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at,
    atransf = exp,
    refline = NA,
    xlab    = "",
    rows    = rows,
    psize   = 1.2,
    cex     = 1.1,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE
  )
  
  usr <- par("usr")
  
  # ---- Zebra shading ----
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(
        usr[1], rows[i] - 0.5,
        usr[2], rows[i] + 0.5,
        col = bg_alt, border = NA
      )
    }
  }
  
  # ---- Header band + line ----
  rect(usr[1], y_header_line - 0.25, usr[2], y_header_line + 0.25,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # ---- Null line at OR=1 (log=0) ----
  segments(0, usr[3], 0, y_header_line, lty = 2, lwd = 1)
  
  # ---- CI lines + caps + arrows when very wide CIs ----
  arrow_dx <- 0.20
  for (i in seq_along(rows)) {
    yi_i  <- df$yi[i]
    sei_i <- sqrt(df$vi[i])
    lcl_i <- yi_i - 1.96 * sei_i
    ucl_i <- yi_i + 1.96 * sei_i
    row_i <- rows[i]
    
    lcl_plot <- max(lcl_i, alim[1])
    ucl_plot <- min(ucl_i, alim[2])
    
    segments(lcl_plot, row_i, ucl_plot, row_i, col = col_ci, lwd = 2.2)
    segments(lcl_plot, row_i - 0.06, lcl_plot, row_i + 0.06, col = col_ci, lwd = 2.2)
    segments(ucl_plot, row_i - 0.06, ucl_plot, row_i + 0.06, col = col_ci, lwd = 2.2)
    
    if (lcl_i < alim[1]) {
      arrows(alim[1], row_i, alim[1] + arrow_dx, row_i,
             length = 0.08, angle = 30, col = col_ci, lwd = 2)
    }
    if (ucl_i > alim[2]) {
      arrows(alim[2], row_i, alim[2] + arrow_dx, row_i,
             length = 0.08, angle = 30, col = col_ci, lwd = 2)
    }
  }
  
  # ---- Squares sized by RE weights ----
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$yi, rows, pch = 15, cex = 1.0 + 1.4 * w_scaled, col = col_sq)
  
  # ---- Column headers ----
  text(x_study,  y_header_line + 0.35, "Study",               font = 2, adj = 0)
  text(x_n,      y_header_line + 0.35, "N",                   font = 2, adj = 1)
  text(x_weight, y_header_line + 0.35, "Weight (%)",          font = 2, adj = 0)
  text(x_est,    y_header_line + 0.35, "OR per 1 SD (95% CI)", font = 2, adj = 0)
  
  # ---- Row labels ----
  text(x_study, rows, df$Study, adj = 0)
  text(x_n, rows, df$N, adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est, rows, fmt_ci(df$OR, df$LCL, df$UCL), adj = 0)
  
  # ---- Random-effects diamond ----
  metafor::addpoly(
    res_re,
    row = row_re,
    atransf = exp,
    mlab = "",
    col = col_dia,
    border = col_dia
  )
  
  # ---- Random-effects label + estimate ----
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- exp(c(res_re$b, res_re$ci.lb, res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "—", adj = 0, cex = 1.10)
  
  # ---- Heterogeneity ----
  text(x_study, row_het, het_txt, adj = 0, cex = 1.10)
}

##################################################################################
# Save in multiple formats
##################################################################################

pdf(file.path(plots_dir, "Forestplot_CS_TMAO_HTN_per1SD.pdf"), width = 14, height = 9)
make_forestplot_crosssectional()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_CS_TMAO_HTN_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot_crosssectional()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_CS_TMAO_HTN_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot_crosssectional()
dev.off()

message("Saved cross-sectional TMAO–HTN forest plots (PDF, TIFF, JPEG) to: ", plots_dir)

################################################################################
### Publication-Ready Diagnostics: Cross-Sectional (TMAO and HTN)
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- Rebuild the cross-sectional HTN-specific dataset and model  ----
df_cs_htn <- read_excel(xlsx_path, sheet = "C-S_TMAO and HTN") %>%
  clean_names() %>%
  mutate(
    OR  = suppressWarnings(as.numeric(dplyr::coalesce(or_per1sd, or))),
    LCL = suppressWarnings(as.numeric(dplyr::coalesce(ci_low_per1sd, or_ci_lower))),
    UCL = suppressWarnings(as.numeric(dplyr::coalesce(ci_high_per1sd, or_ci_upper))),
    N   = suppressWarnings(as.integer(sample_size)),
    Study = paste0(first_author, " (", year, ")"),
    yi  = log(OR),
    sei = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi  = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df_cs_htn) >= 2)

res_re_cs_htn <- metafor::rma(yi, vi, data = df_cs_htn, method = "REML")

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo_cs_htn <- metafor::leave1out(res_re_cs_htn)

if (!"pval" %in% names(res_loo_cs_htn)) {
  res_loo_cs_htn$pval <- 2 * pnorm(abs(res_loo_cs_htn$zval), lower.tail = FALSE)
}

loo_est  <- exp(res_loo_cs_htn$estimate)
loo_lb   <- exp(res_loo_cs_htn$ci.lb)
loo_ub   <- exp(res_loo_cs_htn$ci.ub)
loo_pval <- res_loo_cs_htn$pval
studies  <- df_cs_htn$Study

k <- nrow(df_cs_htn)
y_pos <- rev(seq_len(k))

# Reduced log-scale plotting range for the forest panel
forest_left  <- 0.10
forest_right <- 20

# Columns layout
study_col_x  <- 0.008
effect_col_x <- 35
p_col_x      <- 140
plot_xmax    <- 220

png(file.path(diag_dir, "Sensitivity_LOO_CS_TMAO_HTN.png"),
    width = 15, height = 6.2, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 7.5, 1.2, 5.5), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled OR if study is excluded",
     ylab = "",
     bty = "n",
     log = "x")

axis(1,
     at = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
     labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
     cex.axis = 0.9)

# Study labels
text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

# CIs and points with arrows (very wide CIs)
for (i in seq_len(k)) {
  lb_i  <- loo_lb[i]
  ub_i  <- loo_ub[i]
  est_i <- loo_est[i]
  y_i   <- y_pos[i]
  
  lb_plot <- max(lb_i, forest_left)
  ub_plot <- min(ub_i, forest_right)
  
  segments(lb_plot, y_i, ub_plot, y_i, col = "#19C3A3", lwd = 1.4)
  
  if (lb_i < forest_left) {
    arrows(forest_left, y_i, forest_left * 1.10, y_i,
           length = 0.08, angle = 30, code = 2,
           col = "#19C3A3", lwd = 1.4)
  }
  
  if (ub_i > forest_right) {
    arrows(forest_right / 1.18, y_i, forest_right, y_i,
           length = 0.08, angle = 30, code = 2,
           col = "#19C3A3", lwd = 1.4)
  }
  
  points(min(max(est_i, forest_left), forest_right), y_i,
         pch = 16, cex = 0.9, col = "#19C3A3")
}

# Pooled estimate line
segments(x0 = exp(res_re_cs_htn$b), y0 = 0.5, x1 = exp(res_re_cs_htn$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

# Header rule
segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

# Headers
text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
     adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "OR with 95% CI",
     adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value",
     adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot ----
yi  <- as.numeric(res_re_cs_htn$yi)
sei <- sqrt(res_re_cs_htn$vi)
mu  <- as.numeric(res_re_cs_htn$b)

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
        color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = "Log Odds Ratio",
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_CS_TMAO_HTN.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
egger_cs_htn <- metafor::regtest(res_re_cs_htn, model = "rma", predictor = "sei")
begg_cs_htn  <- metafor::ranktest(res_re_cs_htn)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
cat("\n\n********************************************************\n")
cat("ANALYSIS: CROSS-SECTIONAL - TMAO and HTN (OR)\n")
cat("********************************************************\n")
cat("Number of studies included: ", nrow(df_cs_htn), "\n", sep = "")
cat("Studies included: ", paste(df_cs_htn$Study, collapse = "; "), "\n\n", sep = "")
cat("Egger's Regression Test:\n")
print(egger_cs_htn)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_cs_htn)
sink()

message("Diagnostics for Cross-Sectional TMAO-HTN saved to: ", diag_dir)

##################################################################################
# 3.6. Forest plot for CROSS-SECTIONAL — TMAO and CAD_all (OR)
##################################################################################

#--------------------------------
# Helper functions
#--------------------------------
fmt_ci <- function(or, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", or, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

df_raw <- read_excel(xlsx_path, sheet = "C-S_TMAO and CAD_all") %>%
  clean_names()

df <- df_raw %>%
  mutate(
    OR    = as.numeric(or_per1sd),
    LCL   = as.numeric(ci_low_per1sd),
    UCL   = as.numeric(ci_high_per1sd),
    N     = suppressWarnings(as.integer(sample_size)),
    Study    = paste0(first_author, " (", year, ")"),
    yi    = log(OR),
    sei   = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi    = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

# ---- Random-effects meta-analysis (metafor) ----
res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- meta object (for reporting) ----
m_re <- meta::metagen(
  TE      = df$yi,
  seTE    = sqrt(df$vi),
  studlab = df$Study,
  sm      = "OR",
  comb.fixed  = FALSE,
  comb.random = TRUE,
  method.tau  = "REML",
  hakn        = FALSE
)

# ---- Heterogeneity ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100
tau2 <- as.numeric(res_re$tau2)

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "; I\u00b2=", sprintf("%.1f", I2), "%; Tau\u00b2=", sprintf("%.3f", tau2)
)

# ---- Random-effects weights ----
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

##################################################################################
# Layout 
##################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- Axis limits + ticks ----
x_ticks_major <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 14)

# Minor ticks
x_ticks_minor <- c(seq(0.3, 2.0, by = 0.2),
                   2.5, 3, 4, 5, 6, 8, 10, 12, 14)

x_ticks_minor <- sort(unique(x_ticks_minor[x_ticks_minor >= 0.30 & x_ticks_minor <= 14]))

x_at_major <- log(x_ticks_major)
x_at_minor <- log(x_ticks_minor)

alim <- log(c(0.30, 14))

# ---- Expand forest region ----
xlim <- c(alim[1] - 5.4, alim[2] + 4.0)
ylim <- c(-1.7, k + 2.2)

# ---- Column positions ----
x_study  <- xlim[1] + 0.50
x_n      <- alim[1] - 1.40   # moved left
x_weight <- alim[2] + 0.55
x_est    <- alim[2] + 1.75

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

y_header_line <- k + 1.1

##################################################################################
# Plot function
##################################################################################

make_forestplot_cs <- function() {
  
  par(mar = c(3.8, 6, 2.5, 9))
  par(xpd = NA)
  
  # Base forest: 
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at_major,
    atransf = exp,
    refline = NA,
    xlab    = "",
    rows    = rows,
    psize   = 1.2,
    cex     = 1.1,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE
  )
  
  usr <- par("usr")
  
  # Zebra shading
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(usr[1], rows[i] - 0.5, usr[2], rows[i] + 0.5,
           col = bg_alt, border = NA)
    }
  }
  
  # Header band
  rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # Null line at OR=1 
  segments(0, usr[3], 0, y_header_line, lty = 2)
  
  # Minor ticks (unlabeled)
  axis(
    side   = 1,
    at     = x_at_minor,
    labels = FALSE,
    tcl    = -0.3
  )
  
  # CI lines + caps
  segments(df$yi - 1.96 * sqrt(df$vi), rows,
           df$yi + 1.96 * sqrt(df$vi), rows,
           col = col_ci, lwd = 2.2)
  segments(df$yi - 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi - 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  segments(df$yi + 1.96 * sqrt(df$vi), rows - 0.06,
           df$yi + 1.96 * sqrt(df$vi), rows + 0.06,
           col = col_ci, lwd = 2.2)
  
  # Squares sized by RE weights
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$yi, rows, pch = 15, cex = 1.0 + 1.4 * w_scaled, col = col_sq)
  
  # Column headers
  text(x_study,  k + 1.5, "Study",                font = 2, adj = 0)
  text(x_n,      k + 1.5, "N",                    font = 2, adj = 1)
  text(x_weight, k + 1.5, "Weight (%)",           font = 2, adj = 0)
  text(x_est,    k + 1.5, "OR per 1 SD (95% CI)", font = 2, adj = 0)
  
  # Row labels
  text(x_study, rows, df$Study, adj = 0)
  text(x_n, rows, df$N, adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est, rows, fmt_ci(df$OR, df$LCL, df$UCL), adj = 0)
  
  # Random-effects diamond
  metafor::addpoly(
    res_re,
    row = row_re,
    atransf = exp,
    mlab = "",
    col = col_dia,
    border = col_dia
  )
  
  # Random-effects label + estimate
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- exp(c(res_re$b, res_re$ci.lb, res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "\u2014", adj = 0, cex = 1.10)
  
  # Heterogeneity
  text(x_study, row_het, het_txt, adj = 0, cex = 1.07)
}

##################################################################################
# Save outputs
##################################################################################

pdf(file.path(plots_dir, "Forestplot_CS_TMAO_CAD_per1SD.pdf"), width = 14, height = 9)
make_forestplot_cs()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_CS_TMAO_CAD_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot_cs()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_CS_TMAO_CAD_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot_cs()
dev.off()

message("Saved CROSS-SECTIONAL TMAO–CAD forest plots (PDF, TIFF, JPEG) to: ", plots_dir)

# --- Printing summary ---
print(summary(m_re))

################################################################################
### Publication-Ready Diagnostics: Cross-Sectional (TMAO and CAD)
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- Rebuild the cross-sectional CAD-specific dataset and model  ----
df_cs_cad <- read_excel(xlsx_path, sheet = "C-S_TMAO and CAD_all") %>%
  clean_names() %>%
  mutate(
    OR    = as.numeric(or_per1sd),
    LCL   = as.numeric(ci_low_per1sd),
    UCL   = as.numeric(ci_high_per1sd),
    N     = suppressWarnings(as.integer(sample_size)),
    Study = paste0(first_author, " (", year, ")"),
    yi    = log(OR),
    sei   = (log(UCL) - log(LCL)) / (2 * 1.96),
    vi    = sei^2
  ) %>%
  filter(
    is.finite(OR), is.finite(LCL), is.finite(UCL),
    OR > 0, LCL > 0, UCL > 0,
    is.finite(yi), is.finite(vi), vi > 0
  )

stopifnot(nrow(df_cs_cad) >= 2)

res_re_cs_cad <- metafor::rma(yi, vi, data = df_cs_cad, method = "REML")

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo_cs_cad <- metafor::leave1out(res_re_cs_cad)

if (!"pval" %in% names(res_loo_cs_cad)) {
  res_loo_cs_cad$pval <- 2 * pnorm(abs(res_loo_cs_cad$zval), lower.tail = FALSE)
}

loo_est  <- exp(res_loo_cs_cad$estimate)
loo_lb   <- exp(res_loo_cs_cad$ci.lb)
loo_ub   <- exp(res_loo_cs_cad$ci.ub)
loo_pval <- res_loo_cs_cad$pval
studies  <- df_cs_cad$Study

k <- nrow(df_cs_cad)
y_pos <- rev(seq_len(k))

# Forest region
forest_left  <- min(loo_lb) * 0.96
forest_right <- max(loo_ub) * 1.05

# Columns position
study_col_x  <- forest_left / 1.9
effect_col_x <- forest_right * 1.18
p_col_x      <- forest_right * 1.60
plot_xmax    <- forest_right * 1.90

png(file.path(diag_dir, "Sensitivity_LOO_CS_TMAO_CAD.png"),
    width = 13.5, height = 6.4, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 3.8, 1.2, 3.8), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled OR if study is excluded",
     ylab = "",
     bty = "n",
     log = "x")

axis_ticks <- c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6)
axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)

# Study labels
text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

# CIs and points
for (i in seq_len(k)) {
  lb_i  <- loo_lb[i]
  ub_i  <- loo_ub[i]
  est_i <- loo_est[i]
  y_i   <- y_pos[i]
  
  lb_plot <- max(lb_i, forest_left)
  ub_plot <- min(ub_i, forest_right)
  
  segments(lb_plot, y_i, ub_plot, y_i, col = "#19C3A3", lwd = 1.4)
  
  if (lb_i < forest_left) {
    arrows(forest_left, y_i, forest_left * 1.04, y_i,
           length = 0.08, angle = 30, code = 2,
           col = "#19C3A3", lwd = 1.4)
  }
  
  if (ub_i > forest_right) {
    arrows(forest_right / 1.04, y_i, forest_right, y_i,
           length = 0.08, angle = 30, code = 2,
           col = "#19C3A3", lwd = 1.4)
  }
  
  points(min(max(est_i, forest_left), forest_right), y_i,
         pch = 16, cex = 0.9, col = "#19C3A3")
}

# Pooled estimate line
segments(x0 = exp(res_re_cs_cad$b), y0 = 0.5, x1 = exp(res_re_cs_cad$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

# Header rule
segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

# Headers
text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
     adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "OR with 95% CI",
     adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value",
     adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot ----
yi  <- as.numeric(res_re_cs_cad$yi)
sei <- sqrt(res_re_cs_cad$vi)
mu  <- as.numeric(res_re_cs_cad$b)

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
        color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = "Log Odds Ratio",
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_CS_TMAO_CAD.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
egger_cs_cad <- metafor::regtest(res_re_cs_cad, model = "rma", predictor = "sei")
begg_cs_cad  <- metafor::ranktest(res_re_cs_cad)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
cat("\n\n********************************************************\n")
cat("ANALYSIS: CROSS-SECTIONAL - TMAO and CAD (OR)\n")
cat("********************************************************\n")
cat("Number of studies included: ", nrow(df_cs_cad), "\n", sep = "")
cat("Studies included: ", paste(df_cs_cad$Study, collapse = "; "), "\n\n", sep = "")
cat("Egger's Regression Test:\n")
print(egger_cs_cad)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_cs_cad)
sink()

message("Diagnostics for Cross-Sectional TMAO-CAD saved to: ", diag_dir)

##################################################################################
# 3.8. Forest plot for COHORT studies — Beta coefficients
##################################################################################

#--------------------------------
# Helper functions
#--------------------------------
fmt_ci_beta <- function(beta, lcl, ucl) sprintf("%.3f (%.3f, %.3f)", beta, lcl, ucl)

################################################################################
# Load and prepare data
################################################################################

df_raw <- read_excel(xlsx_path, sheet = "Cohorts_TMAO and cIMT") %>%
  clean_names()

df <- df_raw %>%
  mutate(
    Beta  = as.numeric(beta_per1sd),
    LCL   = as.numeric(ci_low_per1sd),
    UCL   = as.numeric(ci_high_per1sd),
    N     = suppressWarnings(as.integer(sample_size)),
    Study    = paste0(first_author, " (", year, ")"),
    
    # Use CI width to derive SE/variance 
    sei   = (UCL - LCL) / (2 * 1.96),
    sei   = pmax(sei, 1e-6),
    vi    = sei^2,
    
    yi    = Beta
  ) %>%
  filter(
    is.finite(Beta), is.finite(LCL), is.finite(UCL),
    is.finite(vi), vi > 0
  )

stopifnot(nrow(df) >= 2)

################################################################################
# Meta-analysis
################################################################################

res_re <- metafor::rma(yi, vi, data = df, method = "REML")

# ---- Heterogeneity ----
Q    <- res_re$QE
dfQ  <- res_re$k - 1
pQ   <- res_re$QEp
I2   <- max(0, (Q - dfQ) / Q) * 100
tau2 <- as.numeric(res_re$tau2)

het_txt <- paste0(
  "Heterogeneity: Q=", sprintf("%.2f", Q),
  " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ)),
  "; I²=", sprintf("%.1f", I2), "%; Tau²=", sprintf("%.3f", tau2)
)

# ---- Random-effects weights ----
w_re <- 1 / (df$vi + tau2)
df$weight_pct <- round(100 * w_re / sum(w_re), 1)

##################################################################################
# Layout
##################################################################################

k <- nrow(df)
rows <- seq(k, 1)

row_re  <- -0.25
row_het <- -0.85

# ---- Axis limits
range_l <- min(df$LCL, res_re$ci.lb, na.rm = TRUE)
range_u <- max(df$UCL, res_re$ci.ub, na.rm = TRUE)
rng <- range_u - range_l

pad  <- 0.02 * rng 
alim <- c(range_l - pad, range_u + pad)

# ---- Ticks ----
tick_by <- if ((alim[2] - alim[1]) <= 0.7) 0.05 else 0.10

x_ticks_major <- seq(
  from = floor(alim[1] / tick_by) * tick_by,
  to   = ceiling(alim[2] / tick_by) * tick_by,
  by   = tick_by
)
x_at_major <- x_ticks_major

# ---- xlim: modified for fitting text columns ----
xlim <- c(alim[1] - 5.0, alim[2] + 3.2)

# ---- y spacing
row_step <- 1.35
rows <- seq(k, 1) * row_step

row_re  <- -0.25
row_het <- -0.95

ylim <- c(-1.9, max(rows) + 2.2)

# ---- Column positions ----
x_study  <- xlim[1] + 0.45
x_n      <- xlim[1] + 3.60
x_weight <- alim[2] + 0.55
x_est    <- alim[2] + 1.55

# ---- Colors ----
col_ci  <- "#2C73B3"
col_sq  <- "#D55E00"
col_dia <- "#1B9E77"
bg_alt  <- "#F2F4F7"

# ---- Header line position ----
y_header_line <- max(rows) + 1.1

##################################################################################
# Plot function
##################################################################################

make_forestplot_cohorts_cimt_beta <- function() {
  
  par(mar = c(3.8, 6, 2.5, 9))
  par(xpd = NA)
  
  # Base forest 
  metafor::forest(
    res_re,
    slab    = rep("", k),
    xlim    = xlim,
    ylim    = ylim,
    alim    = alim,
    at      = x_at_major,
    refline = NA,
    xlab    = "",
    rows    = rows,
    addfit  = FALSE,
    annot   = FALSE,
    header  = FALSE,
    psize   = 0,
    efac    = 0,
    col     = "white",
    colout  = "white"
  )
  
  usr <- par("usr")
  
  # ---- Zebra shading ----
  for (i in seq_along(rows)) {
    if (i %% 2 == 0) {
      rect(usr[1], rows[i] - 0.5 * row_step, usr[2], rows[i] + 0.5 * row_step,
           col = bg_alt, border = NA)
    }
  }
  
  # ---- Clean header band ----
  rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
       col = "white", border = NA)
  segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
  
  # ---- Null line at Beta = 0 ----
  segments(0, usr[3], 0, y_header_line, lty = 2)
  
  # ---- CI lines + caps + arrows if out of bounds ----
  for (i in seq_along(rows)) {
    row <- rows[i]
    lcl <- df$LCL[i]
    ucl <- df$UCL[i]
    
    lcl_plot <- max(lcl, alim[1])
    ucl_plot <- min(ucl, alim[2])
    
    cap_h     <- 0.12
    arrow_len <- 0.05 * diff(alim)
    
    segments(
      lcl_plot, row, ucl_plot, row,
      col = col_ci, lwd = 2.4, lend = "round"
    )
    
    if (lcl < alim[1]) {
      arrows(
        alim[1], row, alim[1] + arrow_len, row,
        length = 0.10, angle = 25, code = 2,
        col = col_ci, lwd = 2.4
      )
    } else {
      segments(
        lcl_plot, row - cap_h, lcl_plot, row + cap_h,
        col = col_ci, lwd = 2.4, lend = "butt"
      )
    }
    
    if (ucl > alim[2]) {
      arrows(
        alim[2], row, alim[2] - arrow_len, row,
        length = 0.10, angle = 25, code = 2,
        col = col_ci, lwd = 2.4
      )
    } else {
      segments(
        ucl_plot, row - cap_h, ucl_plot, row + cap_h,
        col = col_ci, lwd = 2.4, lend = "butt"
      )
    }
  }
  
  # ---- Squares sized by RE weights ----
  w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
  points(df$Beta, rows, pch = 15, cex = 0.80 + 1.05 * w_scaled, col = col_sq)
  
  # ---- Column headers ----
  text(x_study,  y_header_line + 0.35, "Study",                  font = 2, adj = 0)
  text(x_n,      y_header_line + 0.35, "N",                      font = 2, adj = 1)
  text(x_weight, y_header_line + 0.35, "Weight (%)",             font = 2, adj = 0)
  text(x_est,    y_header_line + 0.35, "Beta per 1 SD (95% CI)", font = 2, adj = 0)
  
  # ---- Row labels ----
  text(x_study,  rows, df$Study, adj = 0)
  text(x_n,      rows, df$N, adj = 1)
  text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0)
  text(x_est,    rows, fmt_ci_beta(df$Beta, df$LCL, df$UCL), adj = 0)
  
  # ---- Random-effects diamond ----
  metafor::addpoly(
    res_re,
    row    = row_re,
    mlab   = "",
    col    = col_dia,
    border = col_dia
  )
  
  # ---- Random-effects label + estimate ----
  text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = 1.10)
  
  re_est <- c(as.numeric(res_re$b), as.numeric(res_re$ci.lb), as.numeric(res_re$ci.ub))
  text(
    x_est, row_re,
    sprintf("%.3f (%.3f, %.3f)", re_est[1], re_est[2], re_est[3]),
    adj = 0, font = 2, cex = 1.10
  )
  text(x_weight, row_re, "—", adj = 0, cex = 1.10)
  
  # ---- Heterogeneity ----
  text(x_study, row_het, het_txt, adj = 0, cex = 1.07)
}

##################################################################################
# Save in multiple formats
##################################################################################

pdf(file.path(plots_dir, "Forestplot_COHORTS_TMAO_cIMT_Beta_per1SD.pdf"), width = 14, height = 9)
make_forestplot_cohorts_cimt_beta()
dev.off()

tiff(
  file.path(plots_dir, "Forestplot_COHORTS_TMAO_cIMT_Beta_per1SD.tiff"),
  width = 14, height = 9, units = "in",
  res = 600, compression = "lzw"
)
make_forestplot_cohorts_cimt_beta()
dev.off()

jpeg(
  file.path(plots_dir, "Forestplot_COHORTS_TMAO_cIMT_Beta_per1SD.jpeg"),
  width = 14, height = 9, units = "in",
  res = 600, quality = 100
)
make_forestplot_cohorts_cimt_beta()
dev.off()

message("Saved COHORT TMAO–cIMT beta forest plots (per 1-SD) to: ", plots_dir)

# --- Printing summary ---
print(summary(res_re))

################################################################################
### Publication-Ready Diagnostics: Cohorts (TMAO and cIMT) - Beta Coefficients
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE)

# ---- Rebuild the cohort cIMT-specific dataset and model  ----
df_cimt <- read_excel(xlsx_path, sheet = "Cohorts_TMAO and cIMT") %>%
  clean_names() %>%
  mutate(
    Beta  = as.numeric(beta_per1sd),
    LCL   = as.numeric(ci_low_per1sd),
    UCL   = as.numeric(ci_high_per1sd),
    N     = suppressWarnings(as.integer(sample_size)),
    Study = paste0(first_author, " (", year, ")"),
    sei   = (UCL - LCL) / (2 * 1.96),
    sei   = pmax(sei, 1e-6),
    vi    = sei^2,
    yi    = Beta
  ) %>%
  filter(
    is.finite(Beta), is.finite(LCL), is.finite(UCL),
    is.finite(vi), vi > 0
  )

stopifnot(nrow(df_cimt) >= 2)

res_re_cimt <- metafor::rma(yi, vi, data = df_cimt, method = "REML")

# ---- 1. Leave-One-Out (Sensitivity) Plot ----
res_loo_cimt <- metafor::leave1out(res_re_cimt)

if (!"pval" %in% names(res_loo_cimt)) {
  res_loo_cimt$pval <- 2 * pnorm(abs(res_loo_cimt$zval), lower.tail = FALSE)
}

loo_est  <- res_loo_cimt$estimate
loo_lb   <- res_loo_cimt$ci.lb
loo_ub   <- res_loo_cimt$ci.ub
loo_pval <- res_loo_cimt$pval
studies  <- df_cimt$Study

k <- nrow(df_cimt)
y_pos <- rev(seq_len(k))

# Linear-scale forest region
forest_left  <- min(loo_lb) - 0.08 * diff(range(c(loo_lb, loo_ub)))
forest_right <- max(loo_ub) + 0.08 * diff(range(c(loo_lb, loo_ub)))
forest_width <- forest_right - forest_left

# Compact layout
study_col_x  <- forest_left - 1.25 * forest_width
effect_col_x <- forest_right + 0.18 * forest_width
p_col_x      <- forest_right + 0.72 * forest_width
plot_xmax    <- forest_right + 1.10 * forest_width

png(file.path(diag_dir, "Sensitivity_LOO_Cohorts_TMAO_cIMT_Beta.png"),
    width = 13.5, height = 6.2, units = "in", res = 600, bg = "white")

par(mar = c(4.2, 4.8, 1.2, 3.6), xpd = NA, bg = "white")

plot(NA, NA,
     xlim = c(study_col_x, plot_xmax),
     ylim = c(0.4, k + 0.95),
     xaxt = "n", yaxt = "n",
     xlab = "Pooled Beta if study is excluded",
     ylab = "",
     bty = "n")

axis_ticks <- pretty(c(forest_left, forest_right))
axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)

# Study labels
text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)

# CI lines and points
segments(loo_lb, y_pos, loo_ub, y_pos, col = "#19C3A3", lwd = 1.4)
points(loo_est, y_pos, pch = 16, cex = 0.9, col = "#19C3A3")

# Pooled estimate line only (red)
segments(x0 = as.numeric(res_re_cimt$b), y0 = 0.5, x1 = as.numeric(res_re_cimt$b), y1 = k + 0.35,
         col = "#E64B5D", lwd = 1)

# Header rule
segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
         lwd = 0.9, col = "black")

# Headers
text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
     adj = 0, font = 2, cex = 0.95)
text(x = effect_col_x, y = k + 0.72, labels = "Beta with 95% CI",
     adj = 0, font = 2, cex = 0.95)
text(x = p_col_x, y = k + 0.72, labels = "p-value",
     adj = 0, font = 2, cex = 0.95)

effect_labels <- sprintf("%.3f [%.3f, %.3f]", loo_est, loo_lb, loo_ub)
p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))

text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)

dev.off()

# ---- 2. Funnel Plot ----
yi  <- as.numeric(res_re_cimt$yi)
sei <- sqrt(res_re_cimt$vi)
mu  <- as.numeric(res_re_cimt$b)

funnel_df <- data.frame(
  effect = yi,
  se = sei
)

se_seq <- seq(0, max(sei) * 1.05, length.out = 500)

funnel_lines <- data.frame(
  se = se_seq,
  left = mu - 1.96 * se_seq,
  right = mu + 1.96 * se_seq
)

p_funnel <- ggplot() +
  geom_line(
    data = funnel_lines,
    aes(x = left, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_line(
    data = funnel_lines,
    aes(x = right, y = se, color = "Pseudo 95% CI"),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
        color = "Estimated \u03b8[IV]"),
    linewidth = 0.7
  ) +
  geom_point(
    data = funnel_df,
    aes(x = effect, y = se, color = "Studies"),
    size = 2.6
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Pseudo 95% CI" = "#D0D0D0",
      "Studies" = "#2C7BE5",
      "Estimated \u03b8[IV]" = "#D81B60"
    ),
    breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
  ) +
  labs(
    title = "Funnel plot",
    x = "Beta coefficient",
    y = "Standard error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12, color = "black")
  )

ggsave(
  filename = file.path(diag_dir, "Funnel_Cohorts_TMAO_cIMT_Beta.png"),
  plot = p_funnel,
  width = 8,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ---- 3. Statistical Tests for Publication Bias ----
egger_cimt <- metafor::regtest(res_re_cimt, model = "rma", predictor = "sei")
begg_cimt  <- metafor::ranktest(res_re_cimt)

sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
cat("\n\n********************************************************\n")
cat("ANALYSIS: COHORTS - TMAO and cIMT (Beta coefficients)\n")
cat("********************************************************\n")
cat("Number of studies included: ", nrow(df_cimt), "\n", sep = "")
cat("Studies included: ", paste(df_cimt$Study, collapse = "; "), "\n\n", sep = "")
cat("Egger's Regression Test:\n")
print(egger_cimt)
cat("\nBegg's Rank Correlation Test:\n")
print(begg_cimt)
sink()

message("Diagnostics for cIMT cohorts saved to: ", diag_dir)

################################################################################
# 3.10. Forest plot for SCFAs-related interventions and SBP, DBP (RCTs)
################################################################################

# ---- Directories (again)
base_dir   <- "~/Documents/PhD files/Thesis/1. SR_GM and CVA/Analysis /Paper1_Metabolites"
output_dir <- file.path(base_dir, "Meta_analysis_outputs")
plots_dir  <- file.path(output_dir, "Forest_plots")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

################################################################################
### Load RCT data
################################################################################
meta_rcts <- readxl::read_excel(xlsx_pathrct) %>% janitor::clean_names()

################################################################################
### Fix publication year
################################################################################
if ("publication_year" %in% names(meta_rcts)) {
  meta_rcts$year_use <- suppressWarnings(as.integer(meta_rcts$publication_year))
} else if ("year" %in% names(meta_rcts)) {
  meta_rcts$year_use <- suppressWarnings(as.integer(meta_rcts$year))
} else {
  meta_rcts$year_use <- NA_integer_
}

################################################################################
### Ensure n(I)/n(C) exists 
################################################################################
meta_rcts <- meta_rcts %>%
  mutate(
    n_interv  = if ("n_interv"  %in% names(.)) suppressWarnings(as.integer(n_interv))  else NA_integer_,
    n_control = if ("n_control" %in% names(.)) suppressWarnings(as.integer(n_control)) else NA_integer_
  ) %>%
  mutate(
    n_interv = ifelse(
      is.na(n_interv),
      dplyr::case_when(
        str_detect(str_to_lower(first_author), "^jama")    ~ 14L,
        str_detect(str_to_lower(first_author), "^li")      ~ 54L,
        str_detect(str_to_lower(first_author), "^verhaar") ~ 11L,
        TRUE ~ NA_integer_
      ),
      n_interv
    ),
    n_control = ifelse(
      is.na(n_control),
      dplyr::case_when(
        str_detect(str_to_lower(first_author), "^jama")    ~ 6L,
        str_detect(str_to_lower(first_author), "^li")      ~ 60L,
        str_detect(str_to_lower(first_author), "^verhaar") ~ 10L,
        TRUE ~ NA_integer_
      ),
      n_control
    )
  )

################################################################################
### Prepare effect-size dataset (WMD in change)
################################################################################
prep_rct_outcome <- function(df, outcome_regex) {
  
  df %>%
    mutate(
      # ---- FIX: prevent duplicated "et al" by cleaning first_author ----
      first_author_clean = str_trim(
        str_replace_all(
          as.character(first_author),
          regex("\\s*,?\\s*et\\s+al\\.?\\s*$", ignore_case = TRUE),
          ""
        )
      ),
      Study = ifelse(
        is.na(year_use),
        paste0(first_author_clean, " et al."),
        paste0(first_author_clean, " et al. (", year_use, ")")
      ),
      CVA = if ("cva_which" %in% names(.)) as.character(cva_which) else NA_character_
    ) %>%
    filter(!is.na(CVA), str_detect(str_to_lower(CVA), str_to_lower(outcome_regex))) %>%
    mutate(
      WMD = suppressWarnings(as.numeric(effect_wmd)),
      SE  = suppressWarnings(as.numeric(se_wmd)),
      LCL = suppressWarnings(as.numeric(ci_lower_wmd)),
      UCL = suppressWarnings(as.numeric(ci_upper_wmd))
    ) %>%
    mutate(
      # If SE missing but CI present, derive SE from CI width
      SE = dplyr::if_else(
        is.na(SE) & is.finite(LCL) & is.finite(UCL),
        (UCL - LCL) / (2 * 1.96),
        SE
      ),
      vi = SE^2,
      yi = WMD
    ) %>%
    mutate(
      N_I = suppressWarnings(as.integer(n_interv)),
      N_C = suppressWarnings(as.integer(n_control))
    ) %>%
    filter(is.finite(yi), is.finite(vi), vi > 0)
}

df_sbp <- prep_rct_outcome(meta_rcts, "sbp")
df_dbp <- prep_rct_outcome(meta_rcts, "dbp")

stopifnot(nrow(df_sbp) >= 2)
stopifnot(nrow(df_dbp) >= 2)

################################################################################
### Meta-analysis (random-effects, REML)
################################################################################
res_sbp <- metafor::rma(yi, vi, data = df_sbp, method = "REML")
res_dbp <- metafor::rma(yi, vi, data = df_dbp, method = "REML")

################################################################################
### Forest plot function
################################################################################
make_rct_forestplot <- function(df, res_re, file_stub) {
  
  # ---- Weights ----
  w_re <- if (!is.null(res_re$wi) && length(res_re$wi) == nrow(df)) res_re$wi else 1 / df$vi
  df$weight_pct <- round(100 * w_re / sum(w_re), 1)
  
  # ---- CI formatter ----
  fmt_ci <- function(est, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", est, lcl, ucl)
  
  # ---- Heterogeneity ----
  Q    <- res_re$QE
  dfQ  <- res_re$k - 1
  pQ   <- res_re$QEp
  I2   <- max(0, (Q - dfQ) / Q) * 100
  tau2 <- res_re$tau2
  
  # ---- Split heterogeneity into 2 lines if text too long ----
  het_line1 <- paste0(
    "Heterogeneity: Q=", sprintf("%.2f", Q),
    " (df=", dfQ, "), p=", ifelse(pQ < 0.001, "<0.001", sprintf("%.3f", pQ))
  )
  het_line2 <- paste0(
    "I2=", sprintf("%.1f", I2), "%; Tau2=", sprintf("%.3f", tau2)
  )
  
  # ---- Layout ----
  k <- nrow(df)
  rows <- seq(k, 1)
  
  row_re  <- -0.25
  row_het <- -0.85
  
  # ---- FIXED EFFECT AXIS LIMITS ----
  alim <- c(-16, 25)
  x_at <- c(-16, -10, -5, 0, 5, 10, 15, 20, 25)
  
  # Extra room for left/right columns
  xlim <- c(alim[1] - 30, alim[2] + 30)
  ylim <- c(-1.7, k + 2.2)
  
  # ---- Column positions ----
  x_study  <- xlim[1] + 1.2
  
  # n(I)/n(C) position 
  x_armn   <- alim[1] - 4.0
  
  # ---- Move weight a bit more to the left ----
  x_weight <- alim[2] + 0.9  
  
  x_est    <- alim[2] + 12.5
  
  # ---- Colors ----
  col_ci  <- "#2C73B3"
  col_sq  <- "#D55E00"
  col_dia <- "#1B9E77"
  bg_alt  <- "#F2F4F7"
  
  y_header_line <- k + 1.1
  
  # ---- Text sizes ----
  cex_txt <- 0.95
  cex_hdr <- 1.00
  cex_het <- 0.95
  
  make_plot <- function() {
    
    par(mar = c(3.8, 12.0, 2.5, 14.0))
    par(xpd = NA)
    
    metafor::forest(
      res_re,
      slab    = rep("", k),
      xlim    = xlim,
      ylim    = ylim,
      alim    = alim,
      at      = x_at,
      refline = 0,
      xlab    = "WMD (mmHg)",
      rows    = rows,
      psize   = 1.2,
      cex     = cex_txt,
      addfit  = FALSE,
      annot   = FALSE,
      header  = FALSE
    )
    
    usr <- par("usr")
    
    # ---- Alternating row shading ----
    for (i in seq_along(rows)) {
      if (i %% 2 == 0) {
        rect(usr[1], rows[i] - 0.5, usr[2], rows[i] + 0.5,
             col = bg_alt, border = NA)
      }
    }
    
    # ---- Clean header band ----
    rect(usr[1], y_header_line - 1.2, usr[2], y_header_line + 0.45,
         col = "white", border = NA)
    segments(usr[1], y_header_line, usr[2], y_header_line, lwd = 1)
    
    # ---- Null line ----
    segments(0, usr[3], 0, y_header_line, lty = 2, lwd = 1)
    
    # ---- CI lines ----
    segments(df$yi - 1.96 * sqrt(df$vi), rows,
             df$yi + 1.96 * sqrt(df$vi), rows,
             col = col_ci, lwd = 2.2)
    segments(df$yi - 1.96 * sqrt(df$vi), rows - 0.06,
             df$yi - 1.96 * sqrt(df$vi), rows + 0.06,
             col = col_ci, lwd = 2.2)
    segments(df$yi + 1.96 * sqrt(df$vi), rows - 0.06,
             df$yi + 1.96 * sqrt(df$vi), rows + 0.06,
             col = col_ci, lwd = 2.2)
    
    # ---- Squares ----
    w_scaled <- sqrt(df$weight_pct / max(df$weight_pct))
    points(df$yi, rows, pch = 15, cex = 0.9 + 1.3 * w_scaled, col = col_sq)
    
    # ---- Column headers ----
    text(x_study,  k + 1.5, "Study",        font = 2, adj = 0, cex = cex_hdr)
    text(x_armn,   k + 1.5, "n(I)/n(C)",    font = 2, adj = 1, cex = cex_hdr)
    text(x_weight, k + 1.5, "Weight (%)",   font = 2, adj = 0, cex = cex_hdr)
    text(x_est,    k + 1.5, "WMD (95% CI)", font = 2, adj = 0, cex = cex_hdr)
    
    # ---- Row labels ----
    text(x_study, rows, df$Study, adj = 0, cex = cex_txt)
    text(x_armn, rows, paste0(df$N_I, " / ", df$N_C), adj = 1, cex = cex_txt)
    text(x_weight, rows, sprintf("%.1f", df$weight_pct), adj = 0, cex = cex_txt)
    text(x_est, rows, fmt_ci(df$yi, df$LCL, df$UCL), adj = 0, cex = cex_txt)
    
    # ---- Random-effects diamond ----
    metafor::addpoly(
      res_re,
      row = row_re,
      mlab = "",
      col = col_dia,
      border = col_dia
    )
    
    # ---- Random-effects label and estimate ----
    text(x_study, row_re, "Random-effects model", adj = 0, font = 2, cex = cex_txt)
    re_est <- c(as.numeric(res_re$b), as.numeric(res_re$ci.lb), as.numeric(res_re$ci.ub))
    text(
      x_est, row_re,
      sprintf("%.2f (%.2f, %.2f)", re_est[1], re_est[2], re_est[3]),
      adj = 0, font = 2, cex = cex_txt
    )
    text(x_weight, row_re, "-", adj = 0, cex = cex_txt)  # ASCII dash
    
    # ---- Heterogeneity (2 lines) ----
    text(x_study, row_het + 0.15, het_line1, adj = 0, cex = cex_het)
    text(x_study, row_het - 0.15, het_line2, adj = 0, cex = cex_het)
  }
  
  # ---- Save outputs ----
  pdf(file.path(plots_dir, paste0(file_stub, ".pdf")), width = 14, height = 9)
  make_plot()
  dev.off()
  
  # macOS quartz: no compression argument
  tiff(
    file.path(plots_dir, paste0(file_stub, ".tiff")),
    width = 14, height = 9, units = "in", res = 600
  )
  make_plot()
  dev.off()
  
  jpeg(
    file.path(plots_dir, paste0(file_stub, ".jpeg")),
    width = 14, height = 9, units = "in", res = 600, quality = 100
  )
  make_plot()
  dev.off()
  
  message("Saved: ", file_stub, " (PDF, TIFF, JPEG) -> ", plots_dir)
}

################################################################################
### Run forest plots (SBP + DBP)
################################################################################
make_rct_forestplot(df_sbp, res_sbp, "Forestplot_RCT_SBP_WMDchange_fixedAxis")
make_rct_forestplot(df_dbp, res_dbp, "Forestplot_RCT_DBP_WMDchange_fixedAxis")

# --- Printing summary ---
print(res_sbp)
print(res_dbp)

################################################################################
### Publication-Ready Diagnostics: RCTs (SBP and DBP)
################################################################################

diag_dir <- file.path(output_dir, "Diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

run_rct_diagnostics <- function(res_re, df, outcome_name, label_prefix) {
  
  # ---- 1. Leave-One-Out (Sensitivity) Plot ----
  res_loo <- metafor::leave1out(res_re)
  
  if (!"pval" %in% names(res_loo)) {
    res_loo$pval <- 2 * pnorm(abs(res_loo$zval), lower.tail = FALSE)
  }
  
  loo_est  <- res_loo$estimate
  loo_lb   <- res_loo$ci.lb
  loo_ub   <- res_loo$ci.ub
  loo_pval <- res_loo$pval
  studies  <- df$Study
  
  k <- nrow(df)
  y_pos <- rev(seq_len(k))
  
  forest_left  <- min(loo_lb) - 0.08 * diff(range(c(loo_lb, loo_ub)))
  forest_right <- max(loo_ub) + 0.08 * diff(range(c(loo_lb, loo_ub)))
  forest_width <- forest_right - forest_left
  
  study_col_x  <- forest_left - 1.35 * forest_width
  effect_col_x <- forest_right + 0.18 * forest_width
  p_col_x      <- forest_right + 0.78 * forest_width
  plot_xmax    <- forest_right + 1.15 * forest_width
  
  file_name_loo <- paste0("Sensitivity_LOO_RCT_", label_prefix, ".png")
  png(file.path(diag_dir, file_name_loo),
      width = 13.5, height = 6.2, units = "in", res = 600, bg = "white")
  
  par(mar = c(4.2, 5.0, 1.2, 3.8), xpd = NA, bg = "white")
  
  plot(NA, NA,
       xlim = c(study_col_x, plot_xmax),
       ylim = c(0.4, k + 0.95),
       xaxt = "n", yaxt = "n",
       xlab = paste("Pooled WMD in", outcome_name, "if study is excluded (mmHg)"),
       ylab = "",
       bty = "n")
  
  axis_ticks <- pretty(c(forest_left, forest_right))
  axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 0.9)
  
  text(x = study_col_x, y = y_pos, labels = studies, adj = 0, cex = 0.92)
  
  segments(loo_lb, y_pos, loo_ub, y_pos, col = "#19C3A3", lwd = 1.4)
  points(loo_est, y_pos, pch = 16, cex = 0.9, col = "#19C3A3")
  
  segments(x0 = as.numeric(res_re$b), y0 = 0.5, x1 = as.numeric(res_re$b), y1 = k + 0.35,
           col = "#E64B5D", lwd = 1)
  
  segments(x0 = study_col_x, y0 = k + 0.5, x1 = plot_xmax, y1 = k + 0.5,
           lwd = 0.9, col = "black")
  
  text(x = study_col_x, y = k + 0.72, labels = "Omitted study",
       adj = 0, font = 2, cex = 0.95)
  text(x = effect_col_x, y = k + 0.72, labels = "WMD with 95% CI",
       adj = 0, font = 2, cex = 0.95)
  text(x = p_col_x, y = k + 0.72, labels = "p-value",
       adj = 0, font = 2, cex = 0.95)
  
  effect_labels <- sprintf("%.2f [%.2f, %.2f]", loo_est, loo_lb, loo_ub)
  p_labels <- ifelse(loo_pval < 0.001, "<0.001", sprintf("%.3f", loo_pval))
  
  text(x = effect_col_x, y = y_pos, labels = effect_labels, adj = 0, cex = 0.9)
  text(x = p_col_x, y = y_pos, labels = p_labels, adj = 0, cex = 0.9)
  
  dev.off()
  
  # ---- 2. Funnel Plot  ----
  funnel_df <- data.frame(
    effect = df$yi,
    se = sqrt(df$vi),
    Study = df$Study
  )
  
  mu <- as.numeric(res_re$b)
  se_seq <- seq(0, max(funnel_df$se) * 1.05, length.out = 500)
  
  funnel_lines <- data.frame(
    se = se_seq,
    left = mu - 1.96 * se_seq,
    right = mu + 1.96 * se_seq
  )
  
  p_funnel <- ggplot() +
    geom_line(
      data = funnel_lines,
      aes(x = left, y = se, color = "Pseudo 95% CI"),
      linewidth = 0.7
    ) +
    geom_line(
      data = funnel_lines,
      aes(x = right, y = se, color = "Pseudo 95% CI"),
      linewidth = 0.7
    ) +
    geom_segment(
      aes(x = mu, xend = mu, y = 0, yend = max(se_seq),
          color = "Estimated \u03b8[IV]"),
      linewidth = 0.7
    ) +
    geom_point(
      data = funnel_df,
      aes(x = effect, y = se, color = "Studies"),
      size = 3
    ) +
    scale_y_reverse(expand = expansion(mult = c(0.01, 0.03))) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Pseudo 95% CI" = "#D0D0D0",
        "Studies" = "#2C7BE5",
        "Estimated \u03b8[IV]" = "#D81B60"
      ),
      breaks = c("Pseudo 95% CI", "Studies", "Estimated \u03b8[IV]")
    ) +
    labs(
      title = "Funnel plot",
      x = paste("Weighted Mean Difference (", outcome_name, ")", sep = ""),
      y = "Standard error"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      panel.grid.major = element_line(color = "#E5E5E5", linetype = "dashed", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 12, color = "black"),
      plot.margin = margin(10, 40, 10, 10)
    )
  
  file_name_funnel <- paste0("Funnel_RCT_", label_prefix, ".png")
  ggsave(
    filename = file.path(diag_dir, file_name_funnel),
    plot = p_funnel,
    width = 9.5,
    height = 6.5,
    dpi = 600,
    bg = "white"
  )
  
  # ---- 3. Statistical Tests ----
  egger_test <- NULL
  begg_test  <- NULL
  
  if (nrow(df) >= 3) {
    egger_test <- tryCatch(
      metafor::regtest(res_re, model = "rma", predictor = "sei"),
      error = function(e) e
    )
    
    begg_test <- tryCatch(
      metafor::ranktest(res_re),
      error = function(e) e
    )
  }
  
  sink(file.path(diag_dir, "Bias_Tests_Results.txt"), append = TRUE)
  cat("\n\n********************************************************\n")
  cat("ANALYSIS: RCT - SCFAs and ", outcome_name, " (WMD)\n", sep = "")
  cat("********************************************************\n")
  cat("Number of studies included: ", nrow(df), "\n", sep = "")
  cat("Studies included: ", paste(df$Study, collapse = "; "), "\n\n", sep = "")
  
  cat("Egger's Regression Test:\n")
  if (nrow(df) < 3) {
    cat("Not performed: fewer than 3 studies were available.\n")
  } else if (inherits(egger_test, "error")) {
    cat("Not performed: ", conditionMessage(egger_test), "\n", sep = "")
  } else {
    print(egger_test)
  }
  
  cat("\nBegg's Rank Correlation Test:\n")
  if (nrow(df) < 3) {
    cat("Not performed: fewer than 3 studies were available.\n")
  } else if (inherits(begg_test, "error")) {
    cat("Not performed: ", conditionMessage(begg_test), "\n", sep = "")
  } else {
    print(begg_test)
  }
  sink()
}

# Run for SBP
run_rct_diagnostics(res_sbp, df_sbp, "SBP", "SBP")

# Run for DBP
run_rct_diagnostics(res_dbp, df_dbp, "DBP", "DBP")

message("Diagnostics for RCTs (SBP & DBP) saved to: ", diag_dir)


############################ End of script #####################################