cascade_model <- function(g, seed_nodes, p = 0.5, nsim = 100, verbose = T) {
  # Independent cascade model 
  # g : network in which simulate
  # seed_nodes : fixed seed vector, or a list with one seed vector per simulation
  # p : probability of infection
  # nsim : number of simulation
  # verbose : parameter to show progress bar
  
  # IDEA:
  # starting from seed_nodes, they are considered ACTIVE, i.e. at the current
  # time step will try to activate, with probability p, its neighbours.
  # Then, each active node become infected.
  
  n <- igraph::vcount(g)
  stopifnot(igraph::is_igraph(g), igraph::is_directed(g), igraph::is_simple(g))
  stopifnot(is.numeric(p), length(p) == 1, is.finite(p), p >= 0, p <= 1,
            is.numeric(nsim), length(nsim) == 1, is.finite(nsim),
            nsim >= 1, nsim == floor(nsim))
  if (is.list(seed_nodes) && length(seed_nodes) != nsim) {
    stop("Provide one seed set for each simulation.")
  }

  validate_seeds <- function(seeds) {
    if (length(seeds) == 0L) return(integer(0))
    if (!is.numeric(seeds) || any(!is.finite(seeds)) ||
        any(seeds != floor(seeds)) || any(seeds < 1 | seeds > n)) {
      stop("Seeds must be vertex indices between 1 and the number of vertices.")
    }
    unique(as.integer(seeds))
  }
  seed_nodes <- if (is.list(seed_nodes)) lapply(seed_nodes, validate_seeds) else validate_seeds(seed_nodes)
  
  # Pre-compute adjacency list once, outside the simulation loop
  adj <- lapply(seq_len(n), function(v) as.integer(igraph::neighbors(g, v, mode = "out")))
  
  result <- vector("list", nsim)
  names(result) <- paste0("iter", seq_len(nsim))
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  Progress [:bar] :percent | ETA: :eta",
      total = nsim, width = 50, clear = FALSE
    )
  }
  
  for (i in seq_len(nsim)) {
    seeds_i <- if (is.list(seed_nodes)) seed_nodes[[i]] else seed_nodes
    
    # Node states: 0 = susceptible, 1 = active, 2 = infected
    # Using a local integer vector is much faster than modifying igraph attributes
    state <- integer(n)
    state[seeds_i] <- 1L
    active <- seeds_i
    t <- 0L
    
    while (length(active) > 0L) {
      new_active <- c()  # nodes activated during this time step
      
      for (v in active) {
        nbrs <- adj[[v]]
        
        # Susceptible neighbors are those with state 0
        susceptible <- nbrs[state[nbrs] == 0L]
        
        if (length(susceptible) > 0L) {
          # Each susceptible neighbor is activated independently with probability p
          activated <- susceptible[runif(length(susceptible)) < p]
          if (length(activated) > 0L) {
            state[activated] <- 1L
            new_active <- c(new_active, activated)
          }
        }
        # v has exhausted its activation attempts: move from active to infected
        state[v] <- 2L
      }
      
      # Only newly activated nodes are active in the next step;
      # duplicates are avoided because we check state == 0 before activating
      active <- new_active
      t <- t + 1L
    }
    
    infected <- which(state == 2L)
    result[[i]] <- list(
      infected   = infected,
      n_infected = length(infected),
      time       = t
    )
    
    if (verbose) pb$tick()
  }
  
  return(list(
    result     = result,
    seed_nodes = seed_nodes,
    nsim       = nsim
  ))
}

