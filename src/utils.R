# Influence in the tree model: the number of edges labeled 1 -> 0.
tree_influence <- function(graph, label) {
  stopifnot(igraph::is_igraph(graph), igraph::is_directed(graph),
            length(label) == igraph::vcount(graph), all(label %in% 0:1))
  edges <- igraph::as_edgelist(graph, names = FALSE)
  sum(label[edges[, 1]] * (1L - label[edges[, 2]]))
}

# Edges use arbitrary character IDs. The nodes file preserves isolated vertices
# and defines the order of the labels returned by the solvers.
read_example_graph <- function(edges_file, nodes_file) {
  edges <- utils::read.csv(edges_file, colClasses = "character")
  nodes <- utils::read.csv(nodes_file, colClasses = "character")
  stopifnot(all(c("from", "to") %in% names(edges)),
            all(c("name", "observed_label") %in% names(nodes)),
            !anyNA(nodes$name), all(nzchar(nodes$name)), !anyDuplicated(nodes$name),
            all(nodes$observed_label %in% c("0", "1")),
            all(edges$from %in% nodes$name), all(edges$to %in% nodes$name))
  nodes$observed_label <- as.integer(nodes$observed_label)
  igraph::graph_from_data_frame(edges[, c("from", "to"), drop = FALSE],
                                directed = TRUE, vertices = nodes)
}

generate_random_directed_tree <- function(n) {
  stopifnot(is.numeric(n), length(n) == 1, is.finite(n), n >= 0, n == floor(n))
  graph <- igraph::make_empty_graph(n, directed = TRUE)
  if (n >= 2) {
    parents <- vapply(2:n, function(child) sample.int(child - 1L, 1L), integer(1))
    graph <- igraph::add_edges(graph, as.vector(rbind(parents, 2:n)))
  }
  graph
}

# Monte Carlo summaries for a fixed seed set. Fractions include the seeds.
summarize_network_influence <- function(x, conf_level = 0.95) {
  stopifnot(is.numeric(x), length(x) >= 1, all(is.finite(x)),
            all(x >= 0 & x <= 1), length(conf_level) == 1,
            is.finite(conf_level), conf_level > 0, conf_level < 1)
  n <- length(x)
  avg <- mean(x)
  std <- stats::sd(x)
  mcse <- std / sqrt(n)
  critical <- if (n > 1) stats::qt((1 + conf_level) / 2, n - 1) else NA_real_
  data.frame(n_sim = n, mean_infected = avg, sd_infected = std, mcse = mcse,
             ci_lower = pmax(0, avg - critical * mcse),
             ci_upper = pmin(1, avg + critical * mcse))
}
