source(file.path("scripts", "setup.R"))
suppressPackageStartupMessages(library(patchwork))

weighted_profiles <- mic %>%
  filter(
    (dataset == "Indicator" & antimicrobial %in% indicator_cluster_panel) |
      (dataset == "3GCR" & antimicrobial %in% gcr_cluster_panel)
  ) %>%
  inner_join(
    memberships,
    by = c("public_isolate_id", "dataset"),
    relationship = "many-to-many"
  ) %>%
  summarise(
    estimated_nwt_probability = weighted.mean(nwt, posterior_probability, na.rm = TRUE),
    .by = c(dataset, cluster, antimicrobial)
  )
write_csv(weighted_profiles, file.path("output", "tables", "cluster_weighted_profiles.csv"))

p_heat <- weighted_profiles %>%
  mutate(
    cluster = factor(cluster, c(cluster_levels$Indicator, cluster_levels$`3GCR`)),
    dataset = factor(dataset, c("Indicator", "3GCR"))
  ) %>%
  ggplot(aes(antimicrobial, cluster, fill = estimated_nwt_probability)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~dataset, scales = "free_y", ncol = 1) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                       labels = percent, name = "estimated\nNWT probability") +
  labs(x = "antimicrobial", y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())

cluster_composition <- isolates %>%
  inner_join(assignments, by = c("public_isolate_id", "dataset")) %>%
  count(dataset, cluster, host, name = "n") %>%
  group_by(dataset, cluster) %>% mutate(proportion = n / sum(n)) %>% ungroup() %>%
  mutate(
    host = factor(host, host_levels),
    dataset = factor(dataset, c("Indicator", "3GCR")),
    cluster = factor(cluster, c(cluster_levels$Indicator, cluster_levels$`3GCR`))
  )
write_csv(cluster_composition, file.path("output", "tables", "cluster_host_composition.csv"))

p_comp <- ggplot(cluster_composition, aes(proportion, cluster, fill = host)) +
  geom_col() + facet_wrap(~dataset, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = host_colours) + scale_x_continuous(labels = percent) +
  labs(x = "host proportion", y = NULL, fill = NULL) + theme_bw(base_size = 8) +
  theme(strip.background = element_blank(), legend.position = "bottom")

save_svg(p_heat | p_comp, "Figure_3_resistance_profile_clusters.svg", 183, 140)

temporal <- isolates %>%
  inner_join(assignments, by = c("public_isolate_id", "dataset")) %>%
  count(dataset, host, year, cluster, name = "n") %>%
  group_by(dataset, host, year) %>%
  filter(sum(n) >= 100) %>% mutate(frequency = n / sum(n)) %>% ungroup() %>%
  mutate(host = factor(host, host_levels))
write_csv(temporal, file.path("output", "tables", "cluster_temporal_frequencies.csv"))

p_temporal <- ggplot(temporal, aes(year, frequency, colour = cluster, group = cluster)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  facet_grid(dataset ~ host) +
  scale_colour_manual(values = cluster_colours) +
  scale_y_continuous(labels = percent) +
  labs(x = "year", y = "within-host cluster frequency", colour = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
        strip.background = element_blank())
save_svg(p_temporal, "Figure_4_temporal_clusters.svg", 183, 120)
