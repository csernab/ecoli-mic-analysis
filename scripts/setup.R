suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

project_dir <- normalizePath(file.path(getwd()), mustWork = TRUE)
if (!file.exists(file.path(project_dir, "data", "ecoli_mic.csv"))) {
  stop("Run scripts from the data_availability directory.")
}

dir.create(file.path(project_dir, "output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "output", "tables"), recursive = TRUE, showWarnings = FALSE)

mic <- read_csv(
  file.path(project_dir, "data", "ecoli_mic.csv"),
  show_col_types = FALSE
)
memberships <- read_csv(
  file.path(project_dir, "data", "cluster_probabilities.csv"),
  show_col_types = FALSE
)
isolates <- mic %>%
  distinct(public_isolate_id, year, host, phenotype_3gcr, dataset)
assignments <- memberships %>%
  filter(assigned) %>%
  transmute(public_isolate_id, dataset, cluster, cluster_max_posterior = posterior_probability)
sales <- read_csv(
  file.path(project_dir, "data", "antimicrobial_sales.csv"),
  show_col_types = FALSE
)

host_levels <- c("Broilers", "Turkeys", "Pigs", "Cattle")
host_colours <- c(
  Broilers = "#E7298A", Turkeys = "#7570B3",
  Pigs = "#D95F02", Cattle = "#1B9E77"
)
dataset_colours <- c(Indicator = "#6387B2", `3GCR` = "#D5C38B")

cluster_levels <- list(
  Indicator = c("I-LowR", "I-FQ", "I-TriR", "I-MDR", "I-MDR+FQ"),
  `3GCR` = c("3G-LowR", "3G-FQ", "3G-MDR", "3G-MDR+FQ")
)
cluster_colours <- c(
  "I-LowR" = "#F4C7A1", "I-FQ" = "#AED581", "I-TriR" = "#F7D6E0",
  "I-MDR" = "#B39DDB", "I-MDR+FQ" = "#7EC8E3",
  "3G-LowR" = "#F4C7A1", "3G-FQ" = "#AED581",
  "3G-MDR" = "#B39DDB", "3G-MDR+FQ" = "#7EC8E3"
)

first_panel <- c("amp", "fot", "taz", "mer", "nal", "cip", "tet",
                 "col", "gen", "tmp", "smx", "chl", "azi", "tgc")
indicator_cluster_panel <- c("amp", "azi", "chl", "cip", "gen", "nal", "smx", "tet", "tmp")
gcr_cluster_panel <- c("azi", "chl", "cip", "gen", "nal", "smx", "tet", "tmp")

save_svg <- function(plot, filename, width, height) {
  ggsave(
    file.path(project_dir, "output", "figures", filename),
    plot, width = width, height = height, units = "mm", dpi = 300,
    device = svglite::svglite
  )
}
