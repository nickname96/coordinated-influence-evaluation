maximize_influence_fixed_k_greedy <- function(tree, k) {
  # Greedy algorithm:
  # 1. Start with all nodes labeled 0.
  # 2. Add k nodes labeled 1, one at a time, always choosing the node
  #    with the largest marginal increase in influence.
  # 3. Repeatedly apply improving switches until no better solution exists.
  #
  # Args:
  #   tree: directed tree as an igraph object
  #   k: number of nodes that must be labeled 1
  #
  # Returns:
  #   A list with:
  #     - label: final 0/1 labeling of the nodes
  #     - influence: total influence of the final 1-nodes
  stopifnot(igraph::is_igraph(tree), is.numeric(k), length(k) == 1,
            is.finite(k), k >= 0, k == floor(k))
  if (!igraph::is_directed(tree) || !igraph::is_simple(tree) ||
      !igraph::is_dag(tree) || any(igraph::degree(tree, mode = "in") > 1)) {
    stop("The input graph must be a forest of rooted directed trees.")
  }
  n <- igraph::vcount(tree)
  if (k > n) {
    message("k is too high. Setting k = number of nodes in the tree.")
    k <- n
  }
  if (k == 0 || k == n || igraph::ecount(tree) == 0) {
    label <- integer(n)
    label[seq_len(k)] <- 1L
    return(list(label = label, influence = 0))
  }
  
  # Initialize all nodes as 0-nodes
  igraph::V(tree)$label <- 0
  
  # Greedily place the first k 1-nodes
  count <- 1
  while (count <= k) {
    dlt <- delta(tree, igraph::V(tree)$label)
    
    # Keep only the nodes with maximum increment
    dlt <- dlt[dlt$increment == max(dlt$increment), ]
    
    # Select the first best candidate
    igraph::V(tree)$label[dlt$node[1]] <- 1
    count <- count + 1
  }
  
  # Repeatedly apply improving switches until convergence
  improved <- TRUE
  while (improved) {
    improved <- FALSE
    
    # Try replacing each current 1-node
    for (j in which(igraph::V(tree)$label == 1)) {
      dlt <- delta(tree, igraph::V(tree)$label, v = j)
      
      if (max(dlt$increment) > 0) {
        improved <- TRUE
        
        # Keep only the nodes with maximum increment
        dlt <- dlt[dlt$increment == max(dlt$increment), ]
        
        # Perform the switch
        igraph::V(tree)$label[j] <- 0
        igraph::V(tree)$label[dlt$node[1]] <- 1
      }
    }
  }
  
  # Compute total influence of the final 1-nodes
  influence <- sum(
    get_degree0(
      tree,
      igraph::V(tree)$label,
      which(igraph::V(tree)$label == 1)
    )
  )
  
  return(list(
    label = igraph::V(tree)$label,
    influence = influence
  ))
}


modified_sample <- function(x) {
  # Return x itself if it has length 1, otherwise sample one element.
  #
  # Args:
  #   x: vector
  #
  # Returns:
  #   A single element from x
  
  if (length(x) == 1) {
    return(x)
  } else {
    return(sample(x, 1))
  }
}


c_l <- function(tree, l, v) {
  # Check whether node v is a child of a 1-node.
  #
  # Args:
  #   tree: directed tree as an igraph object
  #   l: current node labels
  #   v: node index
  #
  # Returns:
  #   1 if v has a parent labeled 1, otherwise 0
  
  if (is.null(v)) {
    return(0)
  }
  
  father <- igraph::neighbors(tree, v, mode = "in")
  
  if (length(father) > 0) {
    return(ifelse(l[father] == 1, 1, 0))
  } else {
    return(0)
  }
}


delta <- function(tree, l, v = NULL, u = which(l == 0)) {
  # Compute the gain in influence obtained by:
  # - adding a 1-node at u, if v is NULL
  # - switching a current 1-node v with a current 0-node u, otherwise
  #
  # Args:
  #   tree: directed tree as an igraph object
  #   l: current node labels
  #   v: node currently labeled 1 to be removed (optional)
  #   u: candidate nodes currently labeled 0
  #
  # Returns:
  #   A data.table with:
  #     - node: candidate node
  #     - increment: gain from adding/switching to that node
  
  dlt <- rep(0, length(u))
  is_child <- c_l(tree, l, v)
  
  index <- 1
  for (i in u) {
    l_new <- l
    
    if (!is.null(v)) {
      l_new[v] <- 0
    }
    
    l_new[i] <- 1
    
    # Degree-0 contribution of the new 1-node i
    d0i <- get_degree0(tree, l_new, i)
    
    # Adjustment:
    # - add 1 if removed node v was itself a child of a 1-node
    # - subtract 1 if added node i becomes a child of a 1-node
    dlt[index] <- d0i + is_child - c_l(tree, l_new, i)
    index <- index + 1
  }
  
  # Contribution of the removed node v; equals 0 when v is NULL
  d0v <- get_degree0(tree, l, v)
  
  return(data.table::data.table(
    node = u,
    increment = dlt - d0v
  ))
}

get_degree0 <- function(tree,l,v = igraph::V(tree)) {
  
  # tree: directed tree (igraph)
  # v : node of tree
  # l : actual label
  
  # return the out-degree on 0-nodes of the vertices in v
  
  if (is.null(v)) {
    
    return(0)
    
  }

  if (length(v) == 0L) return(integer(0))
  
  n1 = which(l == 1) # 1-nodes
  d0 = igraph::degree(tree, v, mode = "out")
  
  for (i in seq_along(d0)) {
    
    neigh = igraph::neighbors(tree,v[i])
    
    # Consider all the neighbours minus the number of 1-nodes
    d0[i] = d0[i] - length(intersect(n1,as.numeric(neigh)))
    
  }
  
  return(d0)
  
}

