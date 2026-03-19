################################################################################
############ Stacked horizontal NOS bar plot (Cross-sectional) #################
################################################################################

# Prepared by: Juliana Hernandez || j.a.hernandezvargas@umcutrecht.nl

if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, ggplot2, forcats, janitor, stringr)

################################################################################
## 1. Paths
################################################################################
base_dir <- "~/Documents/PhD files/Thesis/1. SR_GM and CVA/Analysis /Paper1_Metabolites"
summary_xlsx_path <- file.path(base_dir, "NOS_sectional_input.xlsx")
plots_dir  <- file.path(base_dir, "Meta_analysis_outputs", "Forest_plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

nos_sheet <- "NOS_plot"

################################################################################
## 2. Read + prep
################################################################################
nos_df <- readxl::read_excel(summary_xlsx_path, sheet = nos_sheet) %>%
  janitor::clean_names()

required_cols <- c("first_author", "publication_year", "selection", "comparability", "outcome")
missing_cols <- setdiff(required_cols, names(nos_df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ", paste(missing_cols, collapse = ", "),
    "\nAvailable columns: ", paste(names(nos_df), collapse = ", ")
  )
}

# Compute total_score if missing
if (!("total_score" %in% names(nos_df))) {
  nos_df <- nos_df %>%
    mutate(total_score = selection + comparability + outcome)
}

# Differentiator column so same author/year can appear twice
candidate_id_cols <- c("study", "paper", "title", "journal", "doi", "pmid", "country")
id_col <- candidate_id_cols[candidate_id_cols %in% names(nos_df)][1]

nos_df <- nos_df %>%
  mutate(
    study_label = paste0(first_author, ", et al. (", publication_year, ")"),
    unique_id = if (!is.na(id_col)) {
      paste0(study_label, " | ", .data[[id_col]])
    } else {
      paste0(study_label, " | row_", row_number())
    }
  ) %>%
  arrange(desc(total_score), publication_year, first_author, unique_id) %>%
  mutate(
    # Change bar order (highest to lowest score)
    unique_id = factor(unique_id, levels = rev(unique_id))
  )

################################################################################
## 3. Long format 
################################################################################
nos_long <- nos_df %>%
  pivot_longer(
    cols = c(selection, comparability, outcome),
    names_to = "domain",
    values_to = "stars"
  ) %>%
  mutate(
    domain = recode(domain,
                    selection = "Selection",
                    comparability = "Comparability",
                    outcome = "Outcome"
    ),
    # Legend order
    domain = factor(domain, levels = c("Selection", "Comparability", "Outcome"))
  )

################################################################################
## 4. Plot
################################################################################
p_nos <- ggplot(nos_long, aes(x = stars, y = unique_id, fill = domain)) +
  geom_col(
    width = 0.8,
    position = position_stack(reverse = TRUE)   # <-- KEY FIX
  ) +
  scale_fill_manual(
    breaks = c("Selection", "Comparability", "Outcome"),
    values = c(
      "Selection" = "#4f8fc2",
      "Comparability" = "#f4a259",
      "Outcome" = "#66b266"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 9),
    breaks = 0:9,
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = function(x) str_split_fixed(x, " \\| ", 2)[, 1]) +
  labs(
    x = "NOS stars (0–9)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_blank(),
    legend.position = c(0.86, 0.12),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = "grey80"),
    axis.text.y = element_text(size = 10)
  )

print(p_nos)

################################################################################
## 5. Export
################################################################################
file_base <- file.path(plots_dir, "NOS_CrossSectional_StackedBars_Final")

ggsave(paste0(file_base, ".png"), plot = p_nos, width = 9, height = 10, dpi = 600)
ggsave(paste0(file_base, ".pdf"), plot = p_nos, width = 9, height = 10, device = "pdf")
ggsave(paste0(file_base, ".jpg"), plot = p_nos, width = 14, height = 9, dpi = 600)
ggsave(paste0(file_base, ".tiff"), plot = p_nos, width = 14, height = 9, dpi = 600, compression = "lzw")

message("You've made it: NOS plot saved to: ", plots_dir)

############################ End of script #####################################