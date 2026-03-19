################################################################################
############################ Final Harvest Plot ################################
################################################################################

# Prepared by: Juliana Hernandez || j.a.hernandezvargas@umcutrecht.nl

################################################################################
## 1. Loading libraries
################################################################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, dplyr, stringr, janitor, ggplot2, ggpattern, Cairo)

################################################################################
## 2. Set paths
################################################################################
base_dir <- "~/Documents/PhD files/Thesis/1. SR_GM and CVA/Analysis /Paper1_Metabolites"
summary_xlsx_path <- file.path(base_dir, "Harvest_input_2.0.xlsx")
plots_dir  <- file.path(base_dir, "Meta_analysis_outputs", "Forest_plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

################################################################################
## 3. Prepare Data
################################################################################
plot_ready_df <- readxl::read_excel(summary_xlsx_path) %>%
  janitor::clean_names() %>%
  mutate(
    # ---- Pathway grouping (combine bile acids ----
    pathway_group = case_when(
      str_detect(pathway, regex("bile", ignore_case = TRUE)) ~ "Bile Acids",
      str_detect(pathway, regex("scfa|acetate|propionate|butyrate|valerate|caproic", ignore_case = TRUE)) ~ "SCFAs",
      str_detect(pathway, regex("tmao", ignore_case = TRUE)) ~ "TMAO Precursors",
      str_detect(pathway, regex("tryptophan|indole|kynurenine", ignore_case = TRUE)) ~ "Tryptophan/Indoles",
      TRUE ~ "Others"
    ),
    pathway_group = factor(
      pathway_group,
      levels = c("Bile Acids", "Others", "SCFAs", "TMAO Precursors", "Tryptophan/Indoles")
    ),
    
    # ---- Outcome panels (A=BP/HTN, B=arterial stiffness/PWV, C=atherosclerosis-related)
    panel_label = case_when(
      str_detect(outcome_clean, regex("blood pressure|htn", ignore_case = TRUE)) ~ "A.",
      str_detect(outcome_clean, regex("stiffness|pwv|abi", ignore_case = TRUE)) ~ "B.",
      str_detect(outcome_clean, regex("atherosclerosis|cad|cimt", ignore_case = TRUE)) ~ "C.",
      TRUE ~ outcome_clean
    ),
    panel_label = factor(panel_label, levels = c("A.", "B.", "C.")),
    
    # ---- Standardize RoB (Some concerns -> Moderate) ----
    rob_std = case_when(
      str_to_lower(rob) %in% c("some concerns", "some concern") ~ "Moderate",
      TRUE ~ str_to_title(str_to_lower(rob))
    ),
    rob_std = case_when(
      rob_std %in% c("Low", "Moderate", "High") ~ rob_std,
      TRUE ~ NA_character_
    ),
    
    # ---- Categorise RoB ----
    quality_score = case_when(
      rob_std == "High"     ~ 1,
      rob_std == "Moderate" ~ 2,
      rob_std == "Low"      ~ 3,
      TRUE ~ NA_real_
    ),
    
    # ---- Study design grouping (pattern) ----
    study_design_grouped = case_when(
      str_detect(study_design, regex("nested|case-control", ignore_case = TRUE)) ~ "Case-control / Nested CC",
      str_detect(study_design, regex("cross", ignore_case = TRUE)) ~ "Cross-sectional",
      str_detect(study_design, regex("prospective|cohort|longitud", ignore_case = TRUE)) ~ "Prospective cohort",
      TRUE ~ as.character(study_design)
    ),
    study_design_grouped = factor(
      study_design_grouped,
      levels = c("Case-control / Nested CC", "Cross-sectional", "Prospective cohort")
    ),
    
    # ---- Harmonize association labels // Express everything in terms of direct/inverse (as Oscar suggested) ----
    detailed_assoc_clean = case_when(
      # If dataset uses Direct/Inverse
      str_detect(detailed_assoc, regex("^Significant\\s+Inverse$", ignore_case = TRUE)) ~ "Significant Inverse",
      str_detect(detailed_assoc, regex("^Non-significant\\s+Inverse$", ignore_case = TRUE)) ~ "Non-significant Inverse",
      str_detect(detailed_assoc, regex("^Significant\\s+Direct$", ignore_case = TRUE)) ~ "Significant Direct",
      str_detect(detailed_assoc, regex("^Non-significant\\s+Direct$", ignore_case = TRUE)) ~ "Non-significant Direct",
      
      # If dataset already uses Negative/Positive
      str_detect(detailed_assoc, regex("^Significant\\s+Negative$", ignore_case = TRUE)) ~ "Significant Inverse",
      str_detect(detailed_assoc, regex("^Non-significant\\s+Negative$", ignore_case = TRUE)) ~ "Non-significant Inverse",
      str_detect(detailed_assoc, regex("^Significant\\s+Positive$", ignore_case = TRUE)) ~ "Significant Direct",
      str_detect(detailed_assoc, regex("^Non-significant\\s+Positive$", ignore_case = TRUE)) ~ "Non-significant Direct",
      
      TRUE ~ as.character(detailed_assoc)
    ),
    detailed_assoc_clean = factor(
      detailed_assoc_clean,
      levels = c(
        "Significant Inverse",
        "Non-significant Inverse",
        "Significant Direct",
        "Non-significant Direct"
      )
    ),
    
    # ---- Transparency for significance ----
    is_sig = ifelse(
      str_detect(as.character(detailed_assoc_clean), "^Significant") &
        !str_detect(as.character(detailed_assoc_clean), "Non-significant"),
      1.0, 0.35
    )
  ) %>%
  filter(
    !is.na(panel_label),
    !is.na(detailed_assoc_clean),
    !is.na(quality_score),
    !is.na(pathway_group),
    !is.na(study_design_grouped)
  ) %>%
  group_by(panel_label, detailed_assoc_clean) %>%
  arrange(study_design_grouped, pathway_group, desc(quality_score), first_author, gut_metabolite) %>%
  mutate(x_pos = row_number()) %>%
  ungroup()

################################################################################
## 4. Generate Plot
################################################################################

# Colors by metabolite group
pathway_colors <- c(
  "Bile Acids"         = "#2ca02c",
  "Others"             = "#1f77b4",
  "SCFAs"              = "#9467bd",
  "TMAO Precursors"    = "#d62728",
  "Tryptophan/Indoles" = "#7f7f7f"
)

harvest_p <- ggplot(plot_ready_df, aes(x = x_pos, y = quality_score, fill = pathway_group)) +
  
  geom_col_pattern(
    aes(pattern = study_design_grouped, alpha = is_sig),
    pattern_fill    = "black",
    pattern_color   = "black",
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_size    = 0.3,
    pattern_res     = 600,
    color           = "black",
    size            = 0.2,
    width           = 0.8
  ) +
  
  # Column widths 
  facet_grid(
    panel_label ~ detailed_assoc_clean,
    scales = "free_x",
    space  = "fixed",
    switch = "y"
  ) +
  
  # Colors + patterns
  scale_fill_manual(values = pathway_colors) +
  scale_pattern_manual(values = c("circle", "stripe", "crosshatch")) +
  scale_alpha_identity() +
  
  # RoB axis
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("High", "Mod", "Low"),
    limits = c(0, 4),
    expand = c(0, 0)
  ) +
  
  theme_minimal(base_family = "sans") +
  theme(
    plot.title         = element_blank(),
    plot.subtitle      = element_blank(),
    
    strip.placement    = "outside",
    strip.background   = element_blank(),
    strip.text.y.left  = element_text(face = "bold", size = 14, angle = 0, hjust = 1, margin = margin(r = 15)),
    strip.text.x       = element_text(face = "bold", size = 10, margin = margin(b = 10)),
    
    panel.border       = element_rect(color = "black", fill = NA, size = 0.4),
    panel.spacing.x    = unit(0.4, "lines"),
    panel.spacing.y    = unit(1.5, "lines"),
    
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title         = element_blank(),
    axis.text.y        = element_text(size = 8, color = "black"),
    
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(face = "bold", size = 9),
    legend.text        = element_text(size = 8),
    legend.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(pattern = "none"), order = 1),
    pattern = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(fill = "white", pattern_size = 0.5), order = 2)
  ) +
  labs(fill = "Pathway:", pattern = "Study Design:")

################################################################################
## 5. Export plot
################################################################################

file_base <- file.path(plots_dir, "Harvest_Plot_Final_Refined")

# 1. Save as PNG
ggsave(paste0(file_base, ".png"), plot = harvest_p, width = 14, height = 9, dpi = 600)

# 2. Save as PDF (Standard Device)
ggsave(paste0(file_base, ".pdf"), plot = harvest_p, width = 14, height = 9, device = "pdf")

# 3. Save as JPEG
ggsave(paste0(file_base, ".jpg"), plot = harvest_p, width = 14, height = 9, dpi = 600)

# 4. Save as TIFF
ggsave(paste0(file_base, ".tiff"), plot = harvest_p, width = 14, height = 9, dpi = 600, compression = "lzw")

message("You've made it: Final plots saved.")

############################ End of script #####################################