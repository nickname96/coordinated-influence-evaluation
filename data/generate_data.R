# Run from the repository root. Fixed seeds make the CSV fixtures reproducible.
source("load.R")
dir.create("data", recursive = TRUE, showWarnings = FALSE)
save_graph <- function(graph, prefix, observed) {
  igraph::V(graph)$name <- paste0(prefix, "_", seq_len(igraph::vcount(graph)))
  nodes <- data.frame(name = igraph::V(graph)$name,
                      observed_label = as.integer(seq_len(igraph::vcount(graph)) %in% observed))
  edges <- as.data.frame(igraph::as_edgelist(graph), stringsAsFactors = FALSE)
  names(edges) <- c("from", "to")
  utils::write.csv(nodes, paste0("data/", prefix, "_nodes.csv"), row.names = FALSE)
  utils::write.csv(edges, paste0("data/", prefix, "_edges.csv"), row.names = FALSE)
}

set.seed(42)
tree <- generate_random_directed_tree(24)
save_graph(tree, "tree", c(5, 8, 11, 16, 20))

set.seed(43)
# Three components: two directed trees and one isolated vertex.
forest <- igraph::disjoint_union(generate_random_directed_tree(12),
                                 generate_random_directed_tree(8),
                                 igraph::make_empty_graph(1, directed = TRUE))
save_graph(forest, "forest", c(4, 6, 9, 16))

set.seed(44)
# A directed block network with 3 communities of 30 vertices, plus rings that
# ensure cycles. Within-community links are more likely than cross-community links.
n <- 90L
group <- rep(1:3, each = 30)
probability <- outer(group, group, function(a, b) ifelse(a == b, 0.08, 0.008))
adjacency <- matrix(stats::runif(n * n), n, n) < probability
diag(adjacency) <- FALSE
for (offset in c(0L, 30L, 60L)) {
  from <- offset + seq_len(30)
  to <- offset + c(2:30, 1)
  adjacency[cbind(from, to)] <- TRUE
}
network <- igraph::graph_from_adjacency_matrix(adjacency, mode = "directed")
save_graph(network, "network", c(61, 64, 70, 76, 82))
message("Synthetic tree, forest, and directed network written to data/.")

