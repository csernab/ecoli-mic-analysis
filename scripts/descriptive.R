source(file.path("scripts", "setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(vegan)
})

# Isolate counts and antimicrobial-specific NWT summaries.
isolate_counts <- isolates %>% count(dataset, host, year, name = "n_isolates")
write_csv(isolate_counts, file.path("output", "tables", "isolate_counts.csv"))

nwt_long <- mic %>% select(public_isolate_id, dataset, host, antimicrobial, nwt)

overall_nwt <- nwt_long %>%
  filter(antimicrobial %in% first_panel) %>%
  group_by(dataset, antimicrobial) %>%
  summarise(n = sum(!is.na(nwt)), n_nwt = sum(nwt, na.rm = TRUE),
            proportion_nwt = mean(nwt, na.rm = TRUE), .groups = "drop")
host_nwt <- nwt_long %>%
  filter(antimicrobial %in% first_panel) %>%
  group_by(dataset, host, antimicrobial) %>%
  summarise(n = sum(!is.na(nwt)), n_nwt = sum(nwt, na.rm = TRUE),
            proportion_nwt = mean(nwt, na.rm = TRUE), .groups = "drop")
write_csv(overall_nwt, file.path("output", "tables", "nwt_overall.csv"))
write_csv(host_nwt, file.path("output", "tables", "nwt_by_host.csv"))

p_overall <- ggplot(overall_nwt, aes(antimicrobial, proportion_nwt, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = dataset_colours) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = NULL, y = "NWT proportion", fill = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

selected_ab <- c("amp", "cip", "nal", "tet", "smx", "tmp", "chl", "gen", "azi")
p_host <- host_nwt %>%
  filter(antimicrobial %in% selected_ab) %>%
  mutate(host = factor(host, host_levels)) %>%
  ggplot(aes(antimicrobial, proportion_nwt, colour = host, group = host)) +
  geom_point(size = 1.2) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~dataset, nrow = 1) +
  scale_colour_manual(values = host_colours) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = NULL, y = "NWT proportion", colour = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

save_svg(p_overall / p_host + plot_annotation(tag_levels = "a"),
         "Figure_1_NWT_proportions.svg", 183, 150)

# Resistotype richness using all 14 first-panel WT/NWT calls.
resistotype_data <- mic %>%
  filter(panel == "first") %>%
  select(public_isolate_id, dataset, host, antimicrobial, nwt) %>%
  pivot_wider(names_from = antimicrobial, values_from = nwt) %>%
  unite("resistotype", all_of(first_panel), sep = "", remove = FALSE)
richness <- resistotype_data %>%
  summarise(n_resistotypes = n_distinct(resistotype), .by = dataset)
write_csv(richness, file.path("output", "tables", "resistotype_richness.csv"))

set.seed(19940224)
accum <- list()
for (ds in unique(resistotype_data$dataset)) {
  for (h in host_levels) {
    d <- resistotype_data %>% filter(dataset == ds, host == h)
    incidence <- model.matrix(~ resistotype + 0, data = d)
    a <- specaccum(incidence, method = "random", permutations = 100)
    accum[[paste(ds, h)]] <- tibble(
      dataset = ds, host = h, isolates = a$sites,
      richness = a$richness, sd = a$sd
    )
  }
}
accum <- bind_rows(accum) %>% mutate(host = factor(host, host_levels))
write_csv(accum, file.path("output", "tables", "rarefaction_curves.csv"))

p_rare <- ggplot(accum, aes(isolates, richness, colour = host, fill = host)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.12,
              colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~dataset, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = host_colours) +
  scale_fill_manual(values = host_colours) +
  labs(x = "number of isolates", y = "number of resistotypes", colour = NULL, fill = NULL) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", strip.background = element_blank())
save_svg(p_rare, "Figure_2_resistotype_rarefaction.svg", 183, 100)

# Supplementary descriptive figures.
p_counts <- isolate_counts %>%
  mutate(host = factor(host, host_levels), year = factor(year)) %>%
  ggplot(aes(year, n_isolates, fill = dataset)) +
  geom_col(position = "dodge") + facet_wrap(~host, ncol = 2) +
  scale_fill_manual(values = dataset_colours) +
  labs(x = "year", y = "number of isolates", fill = NULL) +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_svg(p_counts, "Figure_S1_isolate_counts.svg", 150, 105)

p_pheno <- isolates %>% filter(dataset == "3GCR") %>%
  count(host, phenotype_3gcr) %>% group_by(host) %>% mutate(proportion = n / sum(n)) %>%
  ggplot(aes(host, proportion, fill = phenotype_3gcr)) +
  geom_col() + scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "phenotype proportion", fill = NULL) + theme_bw(base_size = 8)
save_svg(p_pheno, "Figure_S2_3GCR_phenotypes.svg", 120, 80)

p_sales <- sales %>%
  mutate(antimicrobial_class = factor(antimicrobial_class)) %>%
  ggplot(aes(year, mg_per_pcu)) + geom_line(linewidth = 0.4) + geom_point(size = 0.7) +
  facet_wrap(~antimicrobial_class, scales = "free_y", ncol = 4) +
  labs(x = "year", y = "mg/PCU") + theme_bw(base_size = 7)
save_svg(p_sales, "Figure_S3_antimicrobial_sales.svg", 183, 125)
