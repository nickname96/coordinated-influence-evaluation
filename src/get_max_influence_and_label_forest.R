get_max_influence_and_label_forest <- function(forest, k, verbose = T) {
  # Compute the optimal influence and labeling for a directed forest
  # with exactly k 1-nodes, allocated jointly across all components.
  stopifnot(igraph::is_igraph(forest))
  stopifnot(is.numeric(k), length(k) == 1, is.finite(k), k >= 0, k == floor(k))
  
  if (!igraph::is_directed(forest) || !igraph::is_simple(forest) || !igraph::is_dag(forest) ||
      any(igraph::degree(forest, mode = "in") > 1)) {
    stop("The input graph must be a forest of rooted directed trees.")
  }
  
  n <- igraph::vcount(forest)
  if (k > n) {
    message("k is too high. Setting k = number of nodes in the forest.")
    k <- n
  }
  
  # These cases have zero influence for every feasible labeling.
  if (k == 0 || k == n || igraph::ecount(forest) == 0) {
    label <- integer(n)
    label[seq_len(k)] <- 1L
    return(list(value = 0, label = label))
  }
  
  roots <- which(igraph::degree(forest, mode = "in") == 0)
  
  # Connect a new dummy root to every component root.
  dummy_root <- n + 1
  tree <- igraph::add_vertices(forest, 1)
  tree <- igraph::add_edges(tree, as.vector(rbind(rep(dummy_root, length(roots)), roots))) 
  dp_res <- influence_constrained_k_binary_tree(tree, k, verbose)
  
  # Fix the dummy root to 0: it uses no 1-node and its edges contribute zero.
  # The DP names vertices by their index. Return only the original vertices,
  # in their original order.
  root_name <- as.character(dummy_root)
  out <- list(
    value = dp_res$S0[root_name, k + 1],
    label = as.integer(dp_res$L0[[root_name]][k + 1, as.character(seq_len(n))])
  )
  
  out
}

