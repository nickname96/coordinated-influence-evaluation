# ============================================================
# Main Function
# ============================================================

influence_constrained_k_binary_tree <- function(tree, k_fix, verbose = FALSE) {
  stopifnot("Input must be an igraph object" = igraph::is_igraph(tree))
  
  n_total <- igraph::vcount(tree)
  if (k_fix > n_total) {
    message("Budget k (", k_fix, ") exceeds tree size. Setting k = ", n_total, ".")
    k_fix <- n_total
  }
  
  igraph::V(tree)$name <- as.character(seq_len(n_total))
  Vnames       <- igraph::V(tree)$name
  
  number_of_nodes_Tv          <- igraph::ego_size(tree, order = n_total, mode = "out")
  names(number_of_nodes_Tv)   <- Vnames
  
  root   <- igraph::V(tree)[number_of_nodes_Tv == n_total]
  leaves <- igraph::V(tree)[number_of_nodes_Tv == 1]
  
  tree_levels <- igraph::distances(tree, v = root, to = igraph::V(tree), mode = "out") |> as.vector()
  depth       <- max(tree_levels)
  
  tensors <- initialize_algorithm(tree, k_fix, n_total, Vnames,
                                  number_of_nodes_Tv, leaves, tree_levels, depth)
  
  S0 <- tensors$S0
  S1 <- tensors$S1
  L0 <- tensors$L0   # now a named list of raw matrices
  L1 <- tensors$L1   # now a named list of raw matrices
  
  if (depth == 1) return(list(S0 = S0, L0 = L0, S1 = S1, L1 = L1))
  
  for (level in seq(depth - 2, 0)) {
    if (verbose) message("Working with level: ", level)
    nodes_lvl <- setdiff(igraph::V(tree)[tree_levels == level], leaves)
    
    for (v in nodes_lvl) {
      v_name   <- Vnames[as.integer(v)]
      v_degree <- igraph::degree(tree, v, mode = "out")
      
      if (v_degree > 2) {
        
        binarized_tree  <- binarize_node(tree, v)
        number_of_nodes_T_binary  <- igraph::ego_size(binarized_tree, order = igraph::vcount(binarized_tree), mode = "out")
        names(number_of_nodes_T_binary) <- igraph::V(binarized_tree)$name
        
        pattern      <- paste0("^", v_name, "_pseudo_[0-9]+$")
        pseudo_nodes <- grep(pattern, igraph::V(binarized_tree)$name, value = TRUE)
        pseudo_nodes <- pseudo_nodes[order(as.integer(stringr::str_extract(pseudo_nodes, "(?<=pseudo_)\\d+")))]
        
        root_b      <- igraph::V(binarized_tree)[number_of_nodes_T_binary == igraph::vcount(binarized_tree)]
        leaves_b    <- igraph::V(binarized_tree)[number_of_nodes_T_binary == 1]
        distances_b <- igraph::distances(binarized_tree, v = root_b, mode = "out")
        names(distances_b) <- igraph::V(binarized_tree)$name
        depth_b <- max(distances_b)
        
        k_fix_b  <- k_fix + v_degree - 2
        tensors_b <- initialize_algorithm(binarized_tree, k_fix_b, igraph::vcount(binarized_tree),
                                          igraph::V(binarized_tree)$name, number_of_nodes_T_binary,
                                          leaves_b, distances_b, max(distances_b))
        
        S0_b <- tensors_b$S0
        S1_b <- tensors_b$S1
        L0_b <- tensors_b$L0
        L1_b <- tensors_b$L1
        
        # Copy already-computed values for nodes common to both trees.
        # Note: v itself was renamed to v_pseudo_1, so it is excluded here.
        nodes_to_update <- intersect(Vnames, igraph::V(binarized_tree)$name)
        k_idx <- 1:(k_fix + 1)
        
        S0_b[nodes_to_update, k_idx] <- S0[nodes_to_update, k_idx]
        S1_b[nodes_to_update, k_idx] <- S1[nodes_to_update, k_idx]
        
        # list-based copy — only copy the columns present in both matrices
        for (vn in nodes_to_update) {
          common_cols <- intersect(colnames(L0[[vn]]), colnames(L0_b[[vn]]))
          if (length(common_cols) > 0) {
            L0_b[[vn]][k_idx, common_cols] <- L0[[vn]][k_idx, common_cols]
            L1_b[[vn]][k_idx, common_cols] <- L1[[vn]][k_idx, common_cols]
          }
        }
        
        n_children  <- 2
        count_start <- 0
        
        for (index in rev(seq_along(pseudo_nodes))) {
          w <- pseudo_nodes[index]
          
          if (distances_b[w] == depth_b - 1) {
            # In this case the node has been already initialized
            count_start <- count_start + 1
            next
          } else if (index == length(pseudo_nodes)) {
            # Base of pseudo-chain: standard binary recurrence
            step_result <- recurrent_step(binarized_tree, k_fix, w,
                                          S0_b, S1_b, L0_b, L1_b,
                                          number_of_nodes_T_binary)
            S0_b <- step_result$S0
            S1_b <- step_result$S1
            L0_b <- step_result$L0
            L1_b <- step_result$L1
            
          } else {
            # Intern of pseudo-chain: apply new recurrence
            children_w       <- igraph::neighbors(binarized_tree, w, mode = "out")
            children_w_names <- igraph::V(binarized_tree)$name[as.integer(children_w)]
            pseudo_child     <- intersect(children_w_names, pseudo_nodes)
            not_pseudo_child <- setdiff(children_w_names, pseudo_nodes)
            
            idx_pseudo <- which(children_w_names == pseudo_child)
            idx_non    <- which(children_w_names == not_pseudo_child)
            
            child_sizes         <- number_of_nodes_T_binary[children_w_names]
            child_sizes_correct <- child_sizes
            child_sizes_correct[idx_pseudo] <- child_sizes_correct[idx_pseudo] - count_start
            
            for (k in seq_len(k_fix)) {
              
              if (k == n_total) {
                S0[v_name, k + 1] <- 0
                S1[v_name, k + 1] <- 0
                L0[[v_name]][k + 1, ] <- as.raw(1L)   # CHANGED
                L1[[v_name]][k + 1, ] <- as.raw(1L)   # CHANGED
                next
              }
              
              # -------- CASE 0: w is a 0-node --------
              if (sum(child_sizes_correct) < k) {
                S0_b[w, k + 1] <- S0_b[w, k]
                L0_b[[w]][k + 1, ] <- L0_b[[w]][k, ]
              } else {
                parts <- rbind(0:k, rev(0:k))
                parts <- parts[, parts[1,] <= child_sizes_correct[1] &
                                 parts[2,] <= child_sizes_correct[2], drop = FALSE]
                
                best_val <- -Inf; best_p <- NULL; best_sel <- NULL
                
                for (i in seq_len(ncol(parts))) {
                  p         <- parts[, i]
                  k_indices <- p + 1
                  val_pseudo  <- S0_b[pseudo_child,     k_indices[idx_pseudo]]
                  S0_val_non  <- S0[not_pseudo_child,   k_indices[idx_non]]
                  S1_val_non  <- S1[not_pseudo_child,   k_indices[idx_non]]
                  sel_non     <- (S1_val_non > S0_val_non) |
                    (p[idx_non] == number_of_nodes_Tv[not_pseudo_child])
                  val_non     <- if (sel_non) S1_val_non else S0_val_non
                  total_val   <- val_pseudo + val_non
                  
                  if (total_val > best_val) {
                    best_val <- total_val; best_p <- p
                    best_sel <- logical(2)
                    best_sel[idx_pseudo] <- FALSE
                    best_sel[idx_non]    <- sel_non
                  }
                }
                
                S0_b[w, k + 1] <- best_val
                
                # named pmax combination of per-node label vectors
                lbl <- as.integer(L0_b[[w]][k + 1, ])
                names(lbl) <- colnames(L0_b[[w]])
                
                pj_p  <- best_p[idx_pseudo] + 1
                src_p <- as.integer(L0_b[[pseudo_child]][pj_p, ])
                names(src_p) <- colnames(L0_b[[pseudo_child]])
                lbl[names(src_p)] <- pmax(lbl[names(src_p)], src_p)
                
                pj_n  <- best_p[idx_non] + 1
                src_n <- as.integer(if (best_sel[idx_non]) L1[[not_pseudo_child]][pj_n, ]
                                    else                   L0[[not_pseudo_child]][pj_n, ])
                names(src_n) <- colnames(L0[[not_pseudo_child]])
                lbl[names(src_n)] <- pmax(lbl[names(src_n)], src_n)
                
                L0_b[[w]][k + 1, ] <- as.raw(lbl)   # CHANGED
              }
              
              # Consolidate CASE 0 back to original tensors when at chain root
              if (index == 1) {
                S0[v_name, min(k_fix, k) + 1] <- S0_b[w, k + 1]
                tmp <- L0_b[[w]][k + 1, ]
                drop_flags <- as.integer(stringr::str_extract(names(tmp), "(?<=_pseudo_)\\d+")) > 1
                tmp <- tmp[!(ifelse(is.na(drop_flags), FALSE, drop_flags))]
                is_p1 <- grepl("_pseudo_1$", names(tmp))
                names(tmp)[is_p1] <- sub("_pseudo_1$", "", names(tmp)[is_p1])
                L0[[v_name]][min(k_fix, k) + 1, names(tmp)] <- tmp # CHANGED
              }
              
              # -------- CASE 1: w is a 1-node --------
              if (sum(child_sizes_correct) == k - 1) {
                # The child must be all 1-nodes since all the pseudo-nodes must be ones
                k_to_select <- k + count_start
                S1_b[w, k_to_select + 1] <- 0
                L1_b[[w]][k_to_select + 1, ] <- as.raw(1L)   # CHANGED
              } else if (sum(child_sizes_correct) < k - 1) {
                # In this case I have less nodes than k to insert
                k_to_select <- k + count_start
                S1_b[w, k_to_select + 1] <- S1_b[w, k_to_select]
                L1_b[[w]][k_to_select + 1, ] <- L1_b[[w]][k_to_select, ]
              } else {
                parts <- rbind(0:(k - 1), rev(0:(k - 1)))
                parts[idx_pseudo, ] <- parts[idx_pseudo, ] + count_start
                parts <- parts[, parts[1,] <= child_sizes[1] &
                                 parts[2,] <= child_sizes[2], drop = FALSE]
                
                best_val <- -Inf; best_p <- NULL; best_sel <- NULL
                k_to_select <- k + count_start
                
                for (i in seq_len(ncol(parts))) {
                  p <- parts[, i]
                  if (any(p > child_sizes)) next
                  k_indices  <- p + 1
                  val_pseudo  <- S1_b[pseudo_child,    k_indices[idx_pseudo]]
                  S0_val_non  <- S0[not_pseudo_child,  k_indices[idx_non]]
                  S1_val_non  <- S1[not_pseudo_child,  k_indices[idx_non]]
                  sel_non     <- (S1_val_non > S0_val_non) |
                    (p[idx_non] == number_of_nodes_Tv[not_pseudo_child])
                  val_non     <- if (sel_non) S1_val_non else S0_val_non
                  total_val   <- val_pseudo + val_non + ifelse(sel_non, 0, 1)
                  
                  if (total_val > best_val) {
                    best_val <- total_val
                    best_p <- p
                    best_sel <- logical(2)
                    best_sel[idx_pseudo] <- TRUE
                    best_sel[idx_non]    <- sel_non
                  }
                }
                
                S1_b[w, k_to_select + 1] <- best_val
                
                # named pmax combination
                lbl <- as.integer(L1_b[[w]][k_to_select + 1, ])
                names(lbl) <- colnames(L1_b[[w]])
                
                pj_p  <- best_p[idx_pseudo] + 1
                src_p <- as.integer(L1_b[[pseudo_child]][pj_p, ])
                names(src_p) <- colnames(L1_b[[pseudo_child]])
                lbl[names(src_p)] <- pmax(lbl[names(src_p)], src_p)
                
                pj_n  <- best_p[idx_non] + 1
                src_n <- as.integer(if (best_sel[idx_non]) L1[[not_pseudo_child]][pj_n, ]
                                    else                   L0[[not_pseudo_child]][pj_n, ])
                names(src_n) <- colnames(L0[[not_pseudo_child]])
                lbl[names(src_n)] <- pmax(lbl[names(src_n)], src_n)
                
                L1_b[[w]][k_to_select + 1, ] <- as.raw(lbl)   # CHANGED
              }
              
              # Consolidate CASE 1 back to original tensors when at chain root
              if (index == 1) {
                S1[v_name, k + 1] <- S1_b[w, k_to_select + 1]
                tmp <- L1_b[[w]][k_to_select + 1, ]
                drop_flags <- as.integer(stringr::str_extract(names(tmp), "(?<=_pseudo_)\\d+")) > 1
                tmp <- tmp[!(ifelse(is.na(drop_flags), FALSE, drop_flags))]
                is_p1 <- grepl("_pseudo_1$", names(tmp))
                names(tmp)[is_p1] <- sub("_pseudo_1$", "", names(tmp)[is_p1])
                L1[[v_name]][k + 1, names(tmp)] <- tmp
              }
            }
          }
          count_start <- count_start + 1
        }
        
      } else {
        step_result <- recurrent_step(tree, k_fix, v, S0, S1, L0, L1, number_of_nodes_Tv)
        S0 <- step_result$S0; S1 <- step_result$S1
        L0 <- step_result$L0; L1 <- step_result$L1
      }
      
      # Free children's L-matrices — they will never be read again.
      # Each child has exactly one parent in a tree; once the parent is processed,
      # L0[[child]] and L1[[child]] serve no further purpose.
      children_of_v <- Vnames[as.integer(igraph::neighbors(tree, v, mode = "out"))]
      for (cn in children_of_v) {
        L0[[cn]] <- NULL
        L1[[cn]] <- NULL
      }
    }
  }
  
  return(list(S0 = S0, L0 = L0, S1 = S1, L1 = L1))
}

# ============================================================
# initialize_algorithm
# ============================================================

initialize_algorithm <- function(tree, k_fix, n_total, Vnames, number_of_nodes_Tv,
                                 leaves, tree_levels, depth) {
  
  S0 <- matrix(0, n_total, k_fix + 1, dimnames = list(Vnames, NULL))
  S1 <- matrix(0, n_total, k_fix + 1, dimnames = list(Vnames, NULL))
  
  # L0/L1 are named lists of raw matrices.
  # L0[[vn]] has dimensions (k_fix+1) x subtree_size(v), columns named by subtree node names.
  # Total memory: O(sum_v subtree_size(v) * k) — O(n log n * k) for balanced trees,
  # vs O(n^2 * k) for the previous monolithic 3-D array.
  subtree_names_list <- setNames(
    lapply(Vnames, function(vn) igraph::subcomponent(tree, vn, mode = "out")$name),
    Vnames
  )
  
  make_L <- function() {
    setNames(lapply(Vnames, function(vn) {
      st <- subtree_names_list[[vn]]
      m  <- matrix(as.raw(0), nrow = k_fix + 1, ncol = length(st))
      colnames(m) <- st
      m
    }), Vnames)
  }
  
  L0 <- make_L()
  L1 <- make_L()
  
  # L1[v, k+1, v] = 1 for k = 1 .. min(subtree_size(v), k_fix)
  for (vn in Vnames) {
    max_k <- min(number_of_nodes_Tv[[vn]], k_fix)
    if (max_k >= 1) L1[[vn]][2:(max_k + 1), vn] <- as.raw(1L)
  }
  
  # Leaves: L0[leaf, 2, leaf] = 1  (mirrors original initialization)
  for (ln in Vnames[as.integer(leaves)]) {
    L0[[ln]][2, ln] <- as.raw(1L)
  }
  
  # Base cases: nodes one level above the leaves
  nodes_lvl <- setdiff(igraph::V(tree)[tree_levels == depth - 1], leaves)
  
  for (v in nodes_lvl) {
    vn            <- Vnames[as.integer(v)]
    children      <- igraph::neighbors(tree, v, mode = "out")
    children_names <- Vnames[as.integer(children)]
    out_deg       <- length(children)
    n_nodes       <- out_deg + 1
    
    # S1 values
    if (n_nodes < k_fix) {
      S1[vn, 2:(n_nodes + 1)] <- seq(out_deg, 0)
    } else {
      S1[vn, 2:(k_fix + 1)] <- seq(out_deg, out_deg - k_fix + 1)
    }
    
    # L1: for k ones, v + first (k-1) children are 1-nodes; v itself already set above
    if (k_fix > 1) {
      for (ki in 2:min(n_nodes, k_fix)) {
        n_incl <- ki - 1
        if (n_incl >= 1 && n_incl <= out_deg)
          L1[[vn]][ki + 1, children_names[1:n_incl]] <- as.raw(1L)
      }
      
      # L0: for k ones, first k children are 1-nodes
      for (ki in 1:min(n_nodes - 1, k_fix)) {
        if (ki <= out_deg)
          L0[[vn]][ki + 1, children_names[1:ki]] <- as.raw(1L)
      }
    }
  }
  return(list(S0 = S0, L0 = L0, S1 = S1, L1 = L1))
}

# ============================================================
# recurrent_step
# ============================================================

recurrent_step <- function(tree, k_fix, v, S0, S1, L0, L1, subtree_sizes) {
  
  # resolve v to both an igraph vertex and a name string (v may arrive as either)
  v_igraph      <- if (is.character(v)) igraph::V(tree)[v] else v
  v_name        <- igraph::V(tree)$name[as.integer(v_igraph)]
  
  children       <- igraph::neighbors(tree, v_igraph, mode = "out")
  n_children     <- length(children)
  children_names <- igraph::V(tree)$name[as.integer(children)]
  
  # name-based subtree_size lookup (safe for pseudo-node names)
  n_nodes     <- subtree_sizes[v_name]
  child_sizes <- subtree_sizes[children_names]
  
  # Pre-compute integer row indices into S0/S1 (avoids repeated string lookup in the hot loop)
  ch_row_idx <- match(children_names, rownames(S0))
  
  for (k in seq_len(min(k_fix, n_nodes))) {
    
    # ======================
    # CASE 0: root is 0-node
    # ======================
    if (k == n_nodes - 1) {
      # derive reachable names from matrix columns instead of igraph::subcomponent()
      reachable_names <- setdiff(colnames(L0[[v_name]]), v_name)
      L0[[v_name]][k + 1, reachable_names] <- as.raw(1L)
      
    } else if (k < n_nodes) {
      
      if (n_children == 2) {
        parts <- rbind(0:k, rev(0:k))
      } else {
        parts <- as.matrix(k)
      }
      
      best_val <- -Inf; best_p <- NULL; best_sel <- NULL
      
      for (i in seq_len(ncol(parts))) {
        p <- parts[, i]
        if (any(p > child_sizes)) next
        idx <- cbind(ch_row_idx, p + 1)
        sel <- (S1[idx] > S0[idx]) | (p == child_sizes)
        val <- sum(ifelse(sel, S1[idx], S0[idx]))
        if (val > best_val) { best_val <- val; best_p <- p; best_sel <- sel }
      }
      
      S0[v_name, k + 1] <- best_val
      
      # read as integer (names preserved from colnames), combine via pmax, store as raw
      lbl <- as.integer(L0[[v_name]][k + 1, ])
      names(lbl) <- colnames(L0[[v_name]])
      
      for (j in seq_along(children_names)) {
        pj  <- best_p[j] + 1
        src <- as.integer(if (best_sel[j]) L1[[children_names[j]]][pj, ]
                          else             L0[[children_names[j]]][pj, ])
        names(src) <- colnames(L1[[children_names[j]]])
        lbl[names(src)] <- pmax(lbl[names(src)], src)
      }
      L0[[v_name]][k + 1, ] <- as.raw(lbl)
    }
    
    # ======================
    # CASE 1: root is 1-node
    # ======================
    if (k == 1) {
      S1[v_name, k + 1] <- igraph::degree(tree, v_igraph, mode = "out")
      next
    } else if (k == n_nodes) {
      S1[v_name, k + 1] <- 0
      L1[[v_name]][k + 1, ] <- as.raw(1L)
      next
    }
    
    # mirror n_children check from CASE 0 (fixes unary-node bug)
    if (n_children == 2) {
      parts <- rbind(0:(k - 1), rev(0:(k - 1)))
    } else {
      parts <- as.matrix(k - 1)
    }
    
    best_val <- -Inf; best_p <- NULL; best_sel <- NULL
    
    for (i in seq_len(ncol(parts))) {
      p <- parts[, i]
      if (any(p > child_sizes)) next
      idx <- cbind(ch_row_idx, p + 1)
      sel <- (S1[idx] > S0[idx]) | (p == child_sizes)
      val <- sum(ifelse(sel, S1[idx], S0[idx])) + sum(!sel)
      if (val > best_val) { best_val <- val; best_p <- p; best_sel <- sel }
    }
    
    S1[v_name, k + 1] <- best_val
    lbl <- as.integer(L1[[v_name]][k + 1, ])
    names(lbl) <- colnames(L1[[v_name]])
    
    for (j in seq_along(children_names)) {
      pj  <- best_p[j] + 1
      src <- as.integer(if (best_sel[j]) L1[[children_names[j]]][pj, ]
                        else             L0[[children_names[j]]][pj, ])
      names(src) <- colnames(L1[[children_names[j]]])
      lbl[names(src)] <- pmax(lbl[names(src)], src)
    }
    L1[[v_name]][k + 1, ] <- as.raw(lbl)
  }
  
  return(list(S0 = S0, L0 = L0, S1 = S1, L1 = L1))
}

# ============================================================
# binarize_node  (unchanged)
# ============================================================

binarize_node <- function(tree, v) {
  if (any(is.null(igraph::V(tree)$name))) igraph::V(tree)$name <- as.character(seq_len(igraph::vcount(tree)))
  igraph::V(tree)$is_pseudo <- FALSE
  if (is.character(v)) v <- igraph::V(tree)[v]
  v <- as.integer(v)
  nodes_to_fix <- v
  while (length(nodes_to_fix) > 0 && any(igraph::degree(tree, nodes_to_fix, mode = "out") > 2)) {
    nodes_to_fix <- nodes_to_fix[igraph::degree(tree, nodes_to_fix, mode = "out") > 2]
    v <- nodes_to_fix[1]
    children  <- igraph::neighbors(tree, v, mode = "out")
    others    <- children[2:length(children)]
    if (!igraph::V(tree)$is_pseudo[v]) igraph::V(tree)$name[v] <- paste0(igraph::V(tree)$name[v], "_pseudo_1")
    v_name    <- igraph::V(tree)$name[v]
    index     <- stringr::str_extract(v_name, "(?<=_pseudo_)[0-9]+$") |> as.integer()
    base_name <- stringr::str_extract(v_name, ".*(?=_pseudo_\\d+$)")
    new_name  <- paste0(base_name, "_pseudo_", index + 1)
    tree      <- igraph::add_vertices(tree, 1, name = new_name)
    new_idx   <- igraph::V(tree)[new_name]
    igraph::V(tree)$is_pseudo[c(v, new_idx)] <- TRUE
    tree <- igraph::delete_edges(tree, igraph::get_edge_ids(tree, as.vector(rbind(v, others))))
    tree <- igraph::add_edges(tree, c(v, new_idx))
    tree <- igraph::add_edges(tree, as.vector(rbind(rep(new_idx, length(others)), others)))
    nodes_to_fix <- c(nodes_to_fix, as.integer(new_idx))
  }
  return(tree)
}
