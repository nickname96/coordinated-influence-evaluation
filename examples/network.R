# Run from the repository root: source("examples/network.R")
source("load.R")
dir.create("results", showWarnings = FALSE)
set.seed(45)

graph <- read_example_graph("data/network_edges.csv", "data/network_nodes.csv")
n <- igraph::vcount(graph)
k <- sum(igraph::V(graph)$observed_label)
nsim <- 100L
probabilities <- seq(0, 1, by = 0.1)

# Sample Random once and keep every seed set fixed across repetitions and
# probabilities. All strategies use the same budget; ranking ties use vertex order.
seed_sets <- list(
  Observed = which(igraph::V(graph)$observed_label == 1L),
  Shell = order(igraph::coreness(graph, mode = "all"), decreasing = TRUE)[seq_len(k)],
  `Out-degree` = order(igraph::degree(graph, mode = "out"), decreasing = TRUE)[seq_len(k)],
  Random = sample.int(n, k)
)
seed_table <- do.call(rbind, lapply(names(seed_sets), function(strategy) {
  data.frame(strategy = strategy, node = igraph::V(graph)$name[seed_sets[[strategy]]])
}))
utils::write.csv(seed_table, "results/network_seeds.csv", row.names = FALSE)

rows <- list()
for (p in probabilities) {
  for (strategy in names(seed_sets)) {
    simulation <- cascade_model(graph, seed_sets[[strategy]], p = p,
                                  nsim = nsim, verbose = FALSE)
    counts <- vapply(simulation$result, `[[`, integer(1), "n_infected")
    rows[[length(rows) + 1L]] <- cbind(data.frame(strategy = strategy, p = p),
                                      summarize_network_influence(counts / n))
  }
}
summary <- do.call(rbind, rows)
utils::write.csv(summary, "results/network_summary.csv", row.names = FALSE)
print(summary[abs(summary$p - 0.3) < 1e-8, ])

summary$strategy <- factor(summary$strategy, levels = c("Random", "Observed", "Shell", "Out-degree"))
network_plot <- ggplot2::ggplot(summary, ggplot2::aes(p, mean_infected, color = strategy, fill = strategy)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = comparison_colors) +
  ggplot2::scale_fill_manual(values = comparison_colors) +
  ggplot2::labs(x = "Infection probability", y = "Mean fraction of activated vertices",
                color = NULL, fill = NULL,
                subtitle = "100 cascades per point; shaded 95% intervals for each fixed seed set") +
  theme_comparison()
save_example_plot(network_plot, "network_influence")
