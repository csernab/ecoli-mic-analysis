source(file.path("scripts", "setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(VGAM)
})

sales_total <- sales %>%
  filter(antimicrobial_class == "Total") %>%
  transmute(year, total_mg_pcu = mg_per_pcu, log2_total = log2(mg_per_pcu))

fit_dataset <- function(dataset_name) {
  levels_cluster <- cluster_levels[[dataset_name]]
  model_data <- isolates %>%
    inner_join(assignments, by = c("public_isolate_id", "dataset")) %>%
    filter(dataset == dataset_name, year <= 2022) %>%
    count(host, year, cluster, name = "n") %>%
    left_join(sales_total, by = "year") %>%
    mutate(
      host = factor(host, host_levels),
      cluster = factor(cluster, levels_cluster)
    )

  fit <- vglm(
    cluster ~ log2_total * host,
    family = multinomial(refLevel = levels_cluster[[1]]),
    weights = n,
    data = model_data
  )

  grid <- expand_grid(
    host = factor(host_levels, host_levels),
    log2_total = seq(min(model_data$log2_total), max(model_data$log2_total), length.out = 100)
  )
  probability <- as_tibble(predict(fit, newdata = grid, type = "response"))
  names(probability) <- levels_cluster
  predictions <- bind_cols(grid, probability) %>%
    pivot_longer(all_of(levels_cluster), names_to = "cluster", values_to = "probability") %>%
    mutate(dataset = dataset_name)

  list(fit = fit, data = model_data, predictions = predictions)
}

extract_sales_or <- function(fit, dataset_name) {
  levels_cluster <- cluster_levels[[dataset_name]]
  comparison_clusters <- levels_cluster[-1]
  coefficients <- coef(fit)
  covariance <- vcov(fit)

  bind_rows(lapply(seq_along(comparison_clusters), function(k) {
    bind_rows(lapply(host_levels, function(h) {
      contrast <- setNames(rep(0, length(coefficients)), names(coefficients))
      base_name <- paste0("log2_total:", k)
      contrast[base_name] <- 1
      if (h != "Broilers") {
        interaction_name <- paste0("log2_total:host", h, ":", k)
        contrast[interaction_name] <- 1
      }
      estimate <- sum(contrast * coefficients)
      se <- sqrt(as.numeric(t(contrast) %*% covariance %*% contrast))
      tibble(
        dataset = dataset_name,
        cluster = comparison_clusters[[k]],
        reference_cluster = levels_cluster[[1]],
        host = h,
        odds_ratio_per_doubling = exp(estimate),
        ci_low = exp(estimate - 1.96 * se),
        ci_high = exp(estimate + 1.96 * se)
      )
    }))
  }))
}

indicator <- fit_dataset("Indicator")
gcr <- fit_dataset("3GCR")
predictions <- bind_rows(indicator$predictions, gcr$predictions)
write_csv(predictions, file.path("output", "tables", "sales_model_predictions.csv"))
write_csv(
  bind_rows(
    extract_sales_or(indicator$fit, "Indicator"),
    extract_sales_or(gcr$fit, "3GCR")
  ),
  file.path("output", "tables", "sales_model_odds_ratios.csv")
)

capture.output(summary(indicator$fit), file = file.path("output", "tables", "sales_model_indicator.txt"))
capture.output(summary(gcr$fit), file = file.path("output", "tables", "sales_model_3GCR.txt"))

plot_predictions <- function(dataset_name, show_legend = TRUE) {
  predictions %>%
    filter(dataset == dataset_name) %>%
    mutate(
      host = factor(host, host_levels),
      cluster = factor(cluster, cluster_levels[[dataset_name]])
    ) %>%
    ggplot(aes(log2_total, probability, colour = host)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~cluster, nrow = 1) +
    scale_colour_manual(values = host_colours) +
    labs(
      x = "log2(total national sales, mg/PCU)",
      y = "predicted probability", colour = NULL,
      title = paste(dataset_name, "E. coli")
    ) +
    theme_bw(base_size = 8) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      strip.background = element_blank(),
      plot.title = element_text(size = 8, face = "bold")
    )
}

p_sales_model <- plot_predictions("Indicator", FALSE) /
  plot_predictions("3GCR", TRUE) +
  plot_layout(heights = c(1, 1), guides = "collect")
save_svg(p_sales_model, "Figure_5_sales_model_probabilities.svg", 183, 115)
