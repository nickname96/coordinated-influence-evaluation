# These settings follow the original research scripts: grey/coral tree nodes,
# theme_bw(), blank facet strips, restrained grids, and legends at the bottom.
comparison_colors <- c("Optimal" = "#d95f02", "Greedy" = "#1b9e77",
                       "Observed" = "#0072B2", "Shell" = "#E69F00",
                       "Out-degree" = "#009E73", "Random" = "#999999")

theme_comparison <- function(base_size = 14) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

plot_labeling <- function(graph, label, title = "Optimal labeling") {
  stopifnot(length(label) == igraph::vcount(graph), all(label %in% 0:1))
  igraph::V(graph)$state <- factor(ifelse(label == 1, "1-node", "0-node"),
                                  levels = c("0-node", "1-node"))
  roots <- which(igraph::degree(graph, mode = "in") == 0)
  ggraph::ggraph(graph, layout = "tree", root = roots) +
    ggraph::geom_edge_link(
      arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
      end_cap = ggraph::circle(3, "mm"), color = "grey70", edge_width = 0.6
    ) +
    ggraph::geom_node_point(ggplot2::aes(color = state), size = 5) +
    ggplot2::scale_color_manual(values = c("0-node" = "grey80", "1-node" = "coral2"),
                                drop = FALSE) +
    ggplot2::labs(title = title,
                  subtitle = paste("Direct influence:", tree_influence(graph, label),
                                   "| Number of 1-nodes:", sum(label)), color = "Node state") +
    ggraph::theme_graph(base_size = 15, base_family = "sans") +
    ggplot2::theme(legend.position = "bottom")
}

save_example_plot <- function(plot, name, width = 8, height = 5) {
  dir.create("figures", showWarnings = FALSE)
  ggplot2::ggsave(file.path("figures", paste0(name, ".pdf")), plot = plot,
                  width = width, height = height, dpi = 300, bg = "white")
  invisible(plot)
}
