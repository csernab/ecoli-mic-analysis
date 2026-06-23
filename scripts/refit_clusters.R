## Optional computationally intensive audit.
## Refit k = 2-10 with 100 deterministic random starts per k directly from
## the released WT/NWT matrix. This script is intentionally not in run_analysis.R.

source(file.path("scripts", "setup.R"))
suppressPackageStartupMessages({
  library(flexmix)
  library(mclust)
  library(purrr)
})

n_restarts <- as.integer(Sys.getenv("N_RESTARTS", "100"))
k_values <- as.integer(strsplit(Sys.getenv("K_VALUES", "2,3,4,5,6,7,8,9,10"), ",")[[1]])
if (is.na(n_restarts) || n_restarts < 2 || anyNA(k_values)) {
  stop("N_RESTARTS must be >=2 and K_VALUES must be comma-separated integers.")
}

mean_pairwise_ari <- function(assignments) {
  if (length(assignments) < 2) return(1)
  pairs <- combn(seq_along(assignments), 2)
  mean(apply(pairs, 2, function(i) {
    adjustedRandIndex(assignments[[i[[1]]]], assignments[[i[[2]]]])
  }))
}

fit_grid <- function(dataset_name, antimicrobials) {
  x <- mic %>%
    filter(dataset == dataset_name) %>%
    filter(antimicrobial %in% antimicrobials) %>%
    select(public_isolate_id, antimicrobial, nwt) %>%
    pivot_wider(names_from = antimicrobial, values_from = nwt) %>%
    arrange(public_isolate_id) %>%
    select(all_of(antimicrobials)) %>%
    as.matrix()

  detailed <- map_dfr(k_values, function(k) {
    message(dataset_name, ": k = ", k)
    map_dfr(seq_len(n_restarts), function(restart) {
      set.seed(restart)
      invisible(capture.output(
        fit <- initFlexmix(
          x ~ 1, k = k, model = FLXMCmvbinary(), nrep = 1,
          control = list(iter.max = 500, minprior = 0.02, verb = 0)
        )
      ))
      post <- posterior(fit)
      max_post <- apply(post, 1, max)
      tibble(
        dataset = dataset_name,
        K = k,
        restart = restart,
        ICL = ICL(fit),
        BIC = BIC(fit),
        mean_max_post = mean(max_post),
        pct_post_ge80 = mean(max_post >= 0.80),
        assignment = list(clusters(fit))
      )
    })
  })

  summary <- detailed %>%
    group_by(dataset, K) %>%
    summarise(
      ICL_best = min(ICL),
      BIC_best = min(BIC),
      mean_max_post_best_ICL = mean_max_post[which.min(ICL)],
      pct_post_ge80_best_ICL = pct_post_ge80[which.min(ICL)],
      ARI_top10_ICL = mean_pairwise_ari(
        assignment[order(ICL)][seq_len(min(10, n()))]
      ),
      .groups = "drop"
    )

  list(detailed = select(detailed, -assignment), summary = summary)
}

indicator <- fit_grid("Indicator", indicator_cluster_panel)
gcr <- fit_grid("3GCR", gcr_cluster_panel)

write_csv(
  bind_rows(indicator$detailed, gcr$detailed),
  file.path("output", "tables", "model_selection_all_refits.csv")
)
write_csv(
  bind_rows(indicator$summary, gcr$summary),
  file.path("output", "tables", "model_selection_refitted_summary.csv")
)
