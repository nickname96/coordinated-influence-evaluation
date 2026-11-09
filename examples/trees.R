# Run from the repository root: source("examples/trees.R")
source("load.R")
dir.create("results", showWarnings = FALSE)

# Solve the included tree and forest with the supplied number of 1-nodes.
labels <- scores <- list()
for (dataset in c("tree", "forest")) {
  graph <- read_example_graph(paste0("data/", dataset, "_edges.csv"),
                              paste0("data/", dataset, "_nodes.csv"))
  observed <- igraph::V(graph)$observed_label
  k <- sum(observed)
  optimal <- get_max_influence_and_label(graph, k, verbose = FALSE)
  greedy <- maximize_influence_fixed_k_greedy(graph, k)
  stopifnot(sum(optimal$label) == k, sum(greedy$label) == k,
            tree_influence(graph, optimal$label) == optimal$value,
            tree_influence(graph, greedy$label) == greedy$influence,
            greedy$influence <= optimal$value)
  labels[[dataset]] <- data.frame(dataset = dataset, node = igraph::V(graph)$name,
                                   observed = observed, optimal = optimal$label,
                                   greedy = greedy$label)
  scores[[dataset]] <- data.frame(dataset = dataset, k = k,
                                   strategy = c("Observed", "Optimal", "Greedy"),
                                   influence = c(tree_influence(graph, observed),
                                                 optimal$value, greedy$influence))
  if (dataset == "tree") {
    tree <- graph
    placement_plot <- plot_labeling(tree, optimal$label)
    save_example_plot(placement_plot, "tree_placement")
  }
}
utils::write.csv(do.call(rbind, labels), "results/tree_labels.csv", row.names = FALSE)
utils::write.csv(do.call(rbind, scores), "results/tree_scores.csv", row.names = FALSE)
print(do.call(rbind, scores))

# Compare the exact and greedy methods on the same tree at several budgets.
budgets <- c(1L, 5L, 10L, 15L, 20L)
comparison <- do.call(rbind, lapply(budgets, function(k) {
  exact_ms <- system.time({
    optimal <- get_max_influence_and_label(tree, k, verbose = FALSE)
  })[["elapsed"]] * 1000
  greedy_ms <- system.time({
    greedy <- maximize_influence_fixed_k_greedy(tree, k)
  })[["elapsed"]] * 1000
  stopifnot(sum(optimal$label) == k, sum(greedy$label) == k,
            tree_influence(tree, optimal$label) == optimal$value,
            tree_influence(tree, greedy$label) == greedy$influence,
            greedy$influence <= optimal$value)
  data.frame(k = k, optimal = optimal$value, greedy = greedy$influence,
               relative_gap = if (optimal$value > 0) (optimal$value - greedy$influence) / optimal$value else NA_real_,
               optimal_ms = exact_ms, greedy_ms = greedy_ms)
}))
utils::write.csv(comparison, "results/tree_comparison.csv", row.names = FALSE)
print(comparison)

plot_data <- rbind(
  data.frame(k = comparison$k, strategy = "Optimal", metric = "Direct influence", value = comparison$optimal),
  data.frame(k = comparison$k, strategy = "Greedy", metric = "Direct influence", value = comparison$greedy),
  data.frame(k = comparison$k, strategy = "Optimal", metric = "Execution time (ms)", value = comparison$optimal_ms),
  data.frame(k = comparison$k, strategy = "Greedy", metric = "Execution time (ms)", value = comparison$greedy_ms)
)
plot_data$strategy <- factor(plot_data$strategy, levels = c("Optimal", "Greedy"))
comparison_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(k, value, color = strategy,
                                                           linetype = strategy, shape = strategy)) +
  ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2.8, stroke = 1) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::scale_color_manual(values = comparison_colors) +
  ggplot2::scale_linetype_manual(values = c("Optimal" = "solid", "Greedy" = "22")) +
  ggplot2::scale_shape_manual(values = c("Optimal" = 16, "Greedy" = 1)) +
  ggplot2::labs(x = "Number of 1-nodes", y = NULL, color = NULL, linetype = NULL, shape = NULL,
                subtitle = "Same synthetic tree and budgets; timings are single-call illustrations") +
  theme_comparison()
save_example_plot(comparison_plot, "exact_vs_greedy", width = 9, height = 4.5)
