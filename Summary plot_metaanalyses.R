################################################################################
################### Summary plot of individual meta-analyses:  #################
################################################################################

# Prepared by: Juliana Hernandez || j.a.hernandezvargas@umcutrecht.nl

################################################################################
### 1. Loading libraries and setting working directory
################################################################################

pacman::p_load(
  readxl, dplyr, stringr, purrr, janitor,
  ggplot2, patchwork, scales, grid
)

# Base directory (Paper 1_Metabolites)
base_dir <- "~/Documents/PhD files/Thesis/1. SR_GM and CVA/Analysis /Paper1_Metabolites"

# Excel input
summary_xlsx_path <- file.path(base_dir, "Summary_meta_analyses_final.xlsx")

# Output directories
output_dir <- file.path(base_dir, "Meta_analysis_outputs")
plots_dir  <- file.path(output_dir, "Forest_plots")

# Create folders if they do not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

################################################################################
### 2. Read summary Excel
################################################################################

sum_df <- readxl::read_excel(summary_xlsx_path) %>%
  janitor::clean_names()

print(names(sum_df))

################################################################################
### 3. Summary forest plot (A/B panels) — pooled meta-analyses table + forest
################################################################################

# Use Helvetica for consistency with other plots
base_family <- "Helvetica"

# Bolder headers
BASE_SIZE <- 11
HDR_FACE  <- "bold"
HDR_COL   <- "black"
HDR_SIZE  <- 3.9 

theme_set(
  theme_void(base_family = base_family) +
    theme(text = element_text(color = "black"))
)

# -------------------------------
# Helper functions
# -------------------------------
is_ratio <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s %in% c("or", "rr", "hr") || grepl("odds|risk|hazard", s)
}

fmt_cases_total <- function(cases, total) {
  if (is.na(total)) return("NA")
  if (is.na(cases)) return(paste0("NA/", formatC(as.integer(total), big.mark = " ", format = "d")))
  paste0(
    formatC(as.integer(cases), big.mark = " ", format = "d"),
    "/",
    formatC(as.integer(total), big.mark = " ", format = "d")
  )
}

fmt_est <- function(effect_size, est, lcl, ucl) {
  if (is_ratio(effect_size)) sprintf("%.2f (%.2f, %.2f)", est, lcl, ucl)
  else sprintf("%.3f (%.3f, %.3f)", est, lcl, ucl)
}

se_from_ci <- function(lcl, ucl) (ucl - lcl) / (2 * 1.96)

compute_pi <- function(effect_size, est, lcl, ucl, tau2) {
  if (any(is.na(c(est, lcl, ucl, tau2)))) return(c(NA_real_, NA_real_))
  se <- se_from_ci(lcl, ucl)
  half <- 1.96 * sqrt(pmax(0, tau2) + se^2)
  c(est - half, est + half)
}

fmt_pi <- function(effect_size, lo, hi) {
  if (is.na(lo) || is.na(hi)) return("NA")
  if (is_ratio(effect_size)) sprintf("%.2f \u2013 %.2f", lo, hi) else sprintf("%.3f \u2013 %.3f", lo, hi)
}

# -------------------------------
# Required columns check
# -------------------------------
stopifnot(all(c(
  "outcome", "exposure_intervention", "n_studies", "study_design",
  "outcome_scale", "pooled_effect", "lower_ci", "upper_ci",
  "i2_percent", "tau2", "effect_size", "grade"
) %in% names(sum_df)))

# -------------------------------
# Robust cases_total creation
# -------------------------------
has_cases_total_col <- "cases_total" %in% names(sum_df)
has_case_num_cols   <- all(c("number_cases", "number_total") %in% names(sum_df))

message(
  "Detected cases/total columns: ",
  paste(
    c(if (has_cases_total_col) "cases_total",
      if (has_case_num_cols) "number_cases+number_total"),
    collapse = ", "
  )
)

if (has_cases_total_col) {
  sum_df$cases_total <- as.character(sum_df$cases_total)
} else if (has_case_num_cols) {
  sum_df$cases_total <- mapply(fmt_cases_total, sum_df$number_cases, sum_df$number_total)
} else {
  sum_df$cases_total <- NA_character_
}

# -------------------------------
# Derived text columns + PI
# -------------------------------
sum_df <- sum_df %>%
  mutate(
    est_txt = mapply(fmt_est, effect_size, pooled_effect, lower_ci, upper_ci),
    i2_txt  = ifelse(is.na(i2_percent), "NA", sprintf("%.1f", i2_percent)),
    t2_txt  = ifelse(is.na(tau2), "NA", sprintf("%.3f", tau2)),
    pi_low  = pmap_dbl(list(effect_size, pooled_effect, lower_ci, upper_ci, tau2),
                       ~ compute_pi(..1, ..2, ..3, ..4, ..5)[1]),
    pi_high = pmap_dbl(list(effect_size, pooled_effect, lower_ci, upper_ci, tau2),
                       ~ compute_pi(..1, ..2, ..3, ..4, ..5)[2]),
    pi_txt  = mapply(fmt_pi, effect_size, pi_low, pi_high)
  )

# Keep Excel order for each panel
dich <- sum_df %>% filter(tolower(outcome_scale) == "dichotomous") %>% mutate(.row = row_number())
cont <- sum_df %>% filter(tolower(outcome_scale) == "continuous")  %>% mutate(.row = row_number())

# -------------------------------
# Plot features
# -------------------------------
COL_CI   <- "#2C73B3"
COL_PT   <- "#D55E00"
ZEBRA    <- "#F2F4F7"
FRAME    <- "#B0B0B0"

# CI 
CI_LW  <- 0.65
CAP_LW <- 0.65
CAP_H  <- 0.10

# Marker size
DIAMOND_SIZE <- 4.7

# Row spacing
ROW_SPACING <- 1.85

# Left columns 
x_out <- 0.05
x_exp <- 0.22
x_n   <- 0.50
x_sd  <- 0.60
x_ct  <- 0.93

# Right columns (enough space between I2/T2/PI)

xr_est <- 0.05
xr_i2  <- 0.42
xr_t2  <- 0.57
xr_pi  <- 0.72
xr_gr  <- 1.05

LEFT_XLIM  <- c(0, 1)
RIGHT_XLIM <- c(0, 1.10)

STD_MARGIN <- margin(6, 6, 6, 6)

# -------------------------------
# Panel builder
# -------------------------------
make_panel <- function(dat, panel_label, xlim, xticks, refline) {
  
  dat <- dat %>%
    arrange(.row) %>%
    mutate(y = rev(seq_len(n())) * ROW_SPACING)
  
  # y-limits fixed for all subplots within a panel
  y_min <- min(dat$y) - ROW_SPACING/2
  y_max <- max(dat$y) + ROW_SPACING*1.20
  
  zebra <- dat %>%
    mutate(is_alt = (.row %% 2 == 0)) %>%
    filter(is_alt) %>%
    transmute(ymin = y - ROW_SPACING/2, ymax = y + ROW_SPACING/2)
  
  dat <- dat %>%
    mutate(
      lcl_plot  = pmax(lower_ci, xlim[1]),
      ucl_plot  = pmin(upper_ci, xlim[2]),
      left_oob  = lower_ci < xlim[1],
      right_oob = upper_ci > xlim[2]
    )
  
  # Header baseline
  y_top <- max(dat$y) + 0.95
  
  # Grey frame limits for the forest panel
  frame_ymin <- y_min
  frame_ymax <- max(dat$y) + ROW_SPACING/2
  
  y_scale_fixed <- scale_y_continuous(
    limits = c(y_min, y_max),
    expand = expansion(mult = 0, add = 0)
  )
  
  # ---- Left table ----
  p_left <- ggplot(dat, aes(y = y)) +
    geom_rect(
      data = zebra,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = ZEBRA, color = NA
    ) +
    annotate("text", x = 0, y = y_top + 0.95, label = panel_label,
             hjust = 0, fontface = "bold", size = 4.4, color = "black") +
    annotate("text", x = x_out, y = y_top, label = "Outcome",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = x_exp, y = y_top, label = "Exposure/Intervention",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = x_n,   y = y_top, label = "n",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 1) +
    annotate("text", x = x_sd,  y = y_top, label = "Study design",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = x_ct,  y = y_top, label = "Cases/Total",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 1) +
    geom_text(aes(x = x_out, label = outcome), hjust = 0, color = "black") +
    geom_text(aes(x = x_exp, label = exposure_intervention), hjust = 0, color = "black") +
    geom_text(aes(x = x_n,   label = n_studies), hjust = 1, color = "black") +
    geom_text(aes(x = x_sd,  label = study_design), hjust = 0, color = "black") +
    geom_text(aes(x = x_ct,  label = cases_total), hjust = 1, color = "black") +
    scale_x_continuous(limits = LEFT_XLIM, expand = c(0, 0)) +
    y_scale_fixed +
    coord_cartesian(clip = "off") +
    theme_void(base_family = base_family, base_size = BASE_SIZE) +
    theme(plot.margin = STD_MARGIN)
  
  # ---- Forest plot ----
  p_forest <- ggplot(dat, aes(y = y)) +
    geom_rect(
      data = zebra,
      aes(xmin = xlim[1], xmax = xlim[2], ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = ZEBRA, color = NA
    ) +
    annotate("rect",
             xmin = xlim[1], xmax = xlim[2],
             ymin = frame_ymin, ymax = frame_ymax,
             fill = NA, color = FRAME, linewidth = 0.7
    ) +
    
    # Null line (manually drawn)
    annotate(
      "segment",
      x = refline, xend = refline,
      y = frame_ymin, yend = frame_ymax,
      linetype = "dashed", linewidth = 0.45, color = "black"
    ) +
    
    geom_segment(aes(x = lcl_plot, xend = ucl_plot, yend = y),
                 color = COL_CI, linewidth = CI_LW) +
    
    geom_segment(
      data = dat %>% filter(!left_oob),
      aes(x = lcl_plot, xend = lcl_plot, y = y - CAP_H, yend = y + CAP_H),
      color = COL_CI, linewidth = CAP_LW
    ) +
    geom_segment(
      data = dat %>% filter(!right_oob),
      aes(x = ucl_plot, xend = ucl_plot, y = y - CAP_H, yend = y + CAP_H),
      color = COL_CI, linewidth = CAP_LW
    ) +
    
    # Arrows for out of bounds CIs
    geom_segment(
      data = dat %>% filter(left_oob),
      aes(x = xlim[1], xend = xlim[1] + 0.03*(xlim[2]-xlim[1]), y = y, yend = y),
      arrow = arrow(type = "closed", length = unit(0.10, "inches")),
      color = COL_CI, linewidth = CI_LW
    ) +
    geom_segment(
      data = dat %>% filter(right_oob),
      aes(x = xlim[2], xend = xlim[2] - 0.03*(xlim[2]-xlim[1]), y = y, yend = y),
      arrow = arrow(type = "closed", length = unit(0.10, "inches")),
      color = COL_CI, linewidth = CI_LW
    ) +
    
    geom_point(aes(x = pooled_effect), shape = 18, size = DIAMOND_SIZE, color = COL_PT) +
    
    scale_x_continuous(limits = xlim, breaks = xticks, expand = c(0, 0)) +
    y_scale_fixed +
    theme_minimal(base_family = base_family, base_size = BASE_SIZE) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black"),
      plot.margin = STD_MARGIN
    )
  
  # ---- Right table ----
  p_right <- ggplot(dat, aes(y = y)) +
    geom_rect(
      data = zebra,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = ZEBRA, color = NA
    ) +
    annotate("text", x = xr_est, y = y_top, label = "Estimate (95% CI)",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = xr_i2,  y = y_top, label = "I\u00b2 (%)",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = xr_t2,  y = y_top, label = "T\u00b2",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = xr_pi,  y = y_top, label = "PI",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 0) +
    annotate("text", x = xr_gr,  y = y_top, label = "GRADE",
             fontface = HDR_FACE, size = HDR_SIZE, color = HDR_COL, hjust = 1) +
    geom_text(aes(x = xr_est, label = est_txt), hjust = 0, color = "black") +
    geom_text(aes(x = xr_i2,  label = i2_txt),  hjust = 0, color = "black") +
    geom_text(aes(x = xr_t2,  label = t2_txt),  hjust = 0, color = "black") +
    geom_text(aes(x = xr_pi,  label = pi_txt),  hjust = 0, color = "black") +
    geom_text(aes(x = xr_gr,  label = grade),   hjust = 1, color = "black") +
    scale_x_continuous(limits = RIGHT_XLIM, expand = c(0, 0)) +
    y_scale_fixed +
    coord_cartesian(clip = "off") +
    theme_void(base_family = base_family, base_size = BASE_SIZE) +
    theme(plot.margin = STD_MARGIN)
  
  # Combine both tables 
  p_left + p_forest + p_right +
    plot_layout(widths = c(5.4, 2.25, 3.9))
}

################################################################################
### 4. Build panels A and B with the correct axes
################################################################################

# Panel A (dichotomous — ratio scale)
panelA <- make_panel(
  dat = dich,
  panel_label = "A.",
  xlim = c(0.5, 3.0),
  xticks = c(0.5, 1, 2, 3),
  refline = 1
)

# Panel B (continuous — linear scale)
mn <- min(cont$lower_ci, na.rm = TRUE)
mx <- max(cont$upper_ci, na.rm = TRUE)
rng <- mx - mn
pad <- ifelse(rng > 0, 0.20 * rng, 1)

xlimB <- c(mn - pad, mx + pad)
xtB   <- pretty(xlimB, n = 6)

panelB <- make_panel(
  dat = cont,
  panel_label = "B.",
  xlim = xlimB,
  xticks = xtB,
  refline = 0
)

final_plot <- panelA / panelB +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(8, 8, 8, 8))

################################################################################
### 5. Save in multiple formats 
################################################################################

out_pdf  <- file.path(plots_dir, "Summary_forestplot_meta_analyses_final_R.pdf")
out_png  <- file.path(plots_dir, "Summary_forestplot_meta_analyses_final_R.png")
out_tiff <- file.path(plots_dir, "Summary_forestplot_meta_analyses_final_R.tiff")
out_jpg  <- file.path(plots_dir, "Summary_forestplot_meta_analyses_final_R.jpeg")

# ---- PDF ----
grDevices::pdf(out_pdf, width = 22, height = 10, onefile = TRUE)
print(final_plot)
grDevices::dev.off()

# ---- PNG ----
ggsave(
  filename = out_png,
  plot     = final_plot,
  width    = 22, height = 10, units = "in",
  dpi      = 300,
  bg       = "white"
)

# ---- TIFF ----
ggsave(
  filename = out_tiff,
  plot     = final_plot,
  width    = 22, height = 10, units = "in",
  dpi      = 600,
  compression = "lzw",
  bg       = "white"
)

# ---- JPEG ----
ggsave(
  filename = out_jpg,
  plot     = final_plot,
  width    = 22, height = 10, units = "in",
  dpi      = 600,
  quality  = 1,
  bg       = "white"
)

message(
  "Saved summary forest plot to:\n",
  out_pdf, "\n",
  out_png, "\n",
  out_tiff, "\n",
  out_jpg
)

############################ End of script #####################################