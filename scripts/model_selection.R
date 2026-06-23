source(file.path("scripts", "setup.R"))

diagnostics <- read_csv(file.path("data", "model_selection.csv"), show_col_types = FALSE)
retained_k <- c(Indicator = 5L, `3GCR` = 4L)

plot_data <- diagnostics %>%
  transmute(
    dataset, K,
    ICL = ICL_best,
    BIC = BIC_best,
    mean_max_post = mean_max_post_best_ICL,
    pct_post_ge80 = pct_post_ge80_best_ICL,
    stability = ARI_top10_ICL,
    retained = K == unname(retained_k[dataset])
  ) %>%
  pivot_longer(c(ICL, BIC, mean_max_post, pct_post_ge80, stability),
               names_to = "metric", values_to = "value") %>%
  mutate(
    dataset = factor(dataset, c("Indicator", "3GCR"), c("Indicator E. coli", "3GCR E. coli")),
    metric = factor(
      metric,
      c("ICL", "BIC", "mean_max_post", "pct_post_ge80", "stability"),
      c("ICL (minimum)", "BIC (minimum)", "Mean maximum posterior",
        "Posterior probability >= 0.80", "Cluster stability (mean ARI)")
    )
  )

p <- ggplot(plot_data, aes(K, value, colour = dataset, group = dataset)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.4) +
  geom_point(data = filter(plot_data, retained), shape = 21, fill = "white",
             size = 3, stroke = 0.9) +
  facet_grid(metric ~ dataset, scales = "free_y") +
  scale_colour_manual(values = c("Indicator E. coli" = "#6387B2", "3GCR E. coli" = "#D5C38B"),
                      guide = "none") +
  scale_x_continuous(breaks = 2:10) +
  labs(x = "number of clusters (k)", y = NULL,
       caption = paste0(
         "ICL and BIC are minima across 100 random restarts per k.\n",
         "Posterior diagnostics use the best ICL-ranked fit; stability is mean ARI among the ten best fits.\n",
         "Open circles mark retained solutions (k = 5 indicator; k = 4 3GCR)."
       )) +
  theme_bw(base_size = 8) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 7, hjust = 0))

save_svg(p, "Figure_S5_model_selection_diagnostics.svg", 165, 190)
