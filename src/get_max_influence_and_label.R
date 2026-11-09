get_max_influence_and_label <- function(tree, k, show_plot = FALSE, verbose = T) {
  # Compute the optimal influence and labeling for a rooted tree or forest
  # with exactly k 1-nodes. Optionally plot the result.
  
  # ---- Dependencies & checks --------------------------------------------
  stopifnot(igraph::is_igraph(tree))
  stopifnot(is.numeric(k), length(k) == 1, is.finite(k), k >= 0, k == floor(k))
  if (!igraph::is_directed(tree) || !igraph::is_simple(tree) ||
      !igraph::is_dag(tree) || any(igraph::degree(tree, mode = "in") > 1)) {
    stop("The input graph must be a forest of rooted directed trees (simple, acyclic, in-degree at most one).")
  }
  if (igraph::vcount(tree) == 0) {
    if (k > 0) message("k is too high. Setting k = number of nodes in the tree.")
    return(list(value = 0, label = integer(0)))
  }
  if(k > igraph::vcount(tree)) {
    message("k is too high. Setting k = number of nodes in the tree.")
    k <- igraph::vcount(tree)
  }
  
  # Identify the roots
  root <- igraph::V(tree)[igraph::degree(tree, mode = "in") == 0]
  is_forest <- !igraph::is_connected(tree, mode = "weak")
  
  # A connected tree must have exactly one root.
  if (!is_forest && length(root) != 1) {
    stop("The input graph is not a rooted tree (unique root not found).")
  }
  
  if (is_forest) {
    out <- get_max_influence_and_label_forest(tree, k, verbose)
  } else if (k == 0 || k == igraph::vcount(tree)) {
    out <- list(value = 0, label = rep(as.integer(k > 0), igraph::vcount(tree)))
  } else {
    # Run dynamic programming
    dp_res <- influence_constrained_k_binary_tree(tree, k, verbose)
    
    # k is the number of 1-nodes, but indices are shifted by +1.
    idx <- k + 1
    root_name <- as.character(as.integer(root))
    node_names <- as.character(seq_len(igraph::vcount(tree)))
    val0 <- dp_res$S0[root_name, idx]
    val1 <- dp_res$S1[root_name, idx]
    
    if (val1 > val0) {
      out <- list(
        value = val1,
        label = as.integer(dp_res$L1[[root_name]][idx, node_names])
      )
    } else {
      out <- list(
        value = val0,
        label = as.integer(dp_res$L0[[root_name]][idx, node_names])
      )
    }
  }
  
  # ---- Optional plot ----------------------------------------------------
  if (isTRUE(show_plot) && igraph::vcount(tree) > 0) {
    
    # Attach node state for plotting (1-node / 0-node)
    igraph::V(tree)$state <- factor(
      ifelse(out$label == 1, "1-node", "0-node"),
      levels = c("0-node", "1-node")
    )
    
    p <- ggraph::ggraph(tree, layout = "tree", root = root) +
      ggraph::geom_edge_link(
        arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
        end_cap = ggraph::circle(3, "mm"),
        color = "grey70", 
        edge_width = 0.6
        ) +
      ggraph::geom_node_point(
        ggplot2::aes(color = state),
        size = 5
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "0-node" = "grey80",
          "1-node" = "coral2"
        )
      ) +
      ggplot2::labs(
        title = paste("Optimal labeling of a", if (is_forest) "forest" else "tree",
                      "with", k, "1-nodes"),
        subtitle = paste("Total influence:", out$value),
        color = "Node state"
      ) +
      ggraph::theme_graph(base_size = 15, base_family = "sans")
    
    print(p)
  }
  
  out
}

