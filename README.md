This repository contains the R code for the algorithms presented in [**Optimal and heuristic strategies for evaluating the influence of coordinated behavior in information cascades and retweet networks**](https://arxiv.org/abs/2609.24398).

## Run the examples

Open R in this directory. Install the dependencies once:

``` r
install.packages(c("igraph", "stringr", "data.table", "progress", "ggplot2", "ggraph"))
```

Then run either example:

``` r
source("examples/trees.R")
source("examples/network.R")
```

- **Trees:** finds exact and greedy placements on a synthetic tree and forest, and compares their influence and execution time at several budgets on the tree.
- **Network:** runs independent cascades on a synthetic directed network, comparing reference, shell, out-degree, and random seed placements.

The scripts create CSV tables in `results/` and PDF plots in `figures/`.

## Code and data

- `src/`: exact dynamic program, forest wrapper, greedy method, cascade simulator, and shared reading and plotting functions.
- `data/`: edge and node CSV files for a 24-node tree, a 21-node forest, and a 90-node directed network. All identifiers and reference labels are synthetic.
- `data/generate_data.R`: recreates the six data files using fixed random seeds.

Edges use `from,to` columns and point in the direction of propagation. Node files use `name,observed_label` columns; labels are 0 or 1. The node table includes isolated vertices and defines the order of returned labels.

### Functions in `src/`

The main entry points are:

| Function | Purpose |
|------------------------------------|------------------------------------|
| `get_max_influence_and_label()` | Finds the exact maximum influence and a corresponding binary labeling for a tree or forest; can also plot it. |
| `get_max_influence_and_label_forest()` | Finds the exact solution for a forest, allocating one total budget across its components. |
| `maximize_influence_fixed_k_greedy()` | Builds a placement greedily, then improves it through swaps; returns its labels and influence. |
| `cascade_model()` | Simulates independent cascades from supplied seed nodes and returns activated vertices and counts for each repetition. |

Utilities in `utils.R` and plotting functions in `plots.R`:

| Function | Purpose |
|------------------------------------|------------------------------------|
| `tree_influence()` | Counts edges from 1-nodes to 0-nodes for a supplied labeling. |
| `read_example_graph()` | Reads edge and node CSV files into a directed graph, preserving node order and reference labels. |
| `generate_random_directed_tree()` | Generates a random tree by attaching each new node to a randomly chosen earlier node. |
| `summarize_network_influence()` | Summarizes activated fractions with their mean, standard deviation, Monte Carlo standard error, and confidence interval. |
| `theme_comparison()` | Applies the shared theme used in the comparison plots. |
| `plot_labeling()` | Draws a tree or forest with nodes colored by their binary labels. |
| `save_example_plot()` | Saves an example plot as a PDF in `figures/`. |

The algorithm files also contain these helpers that are used by the entry points above rather than be called directly:

| Function | Purpose |
|------------------------------------|------------------------------------|
| `influence_constrained_k_binary_tree()` | Runs the exact dynamic program, using auxiliary binary branches for nodes with many children. |
| `initialize_algorithm()` | Initializes the dynamic programming score tables, labeling storage, and base cases. |
| `make_L()` | Local helper inside `initialize_algorithm()` that allocates labeling matrices for each subtree. |
| `recurrent_step()` | Updates a node's dynamic programming tables from the results for its one or two children. |
| `binarize_node()` | Replaces a high-degree branching point with a binary chain of auxiliary nodes for the dynamic program. |
| `get_degree0()` | Counts the outgoing neighbors labeled 0 for each requested node. |
| `c_l()` | Checks whether a node has a parent labeled 1. |
| `delta()` | Computes the change in influence from adding a selected node or swapping a selected and unselected node. |
| `modified_sample()` | Returns a sole element unchanged, or randomly chooses one element from a longer vector. |
| `validate_seeds()` | Local helper inside `cascade_model()` that checks seed indices and removes duplicates. |

## Use your own data

The algorithms use directed `igraph` objects. To load your data, prepare **two CSV files for each graph**:

**Edge file** — each row is a directed edge from `from` to `to`:

``` csv
from,to
A,B
A,C
```

**Node file** — list every node once, including isolated nodes:

``` csv
name,observed_label
A,1
B,0
C,0
```

Some simple observations to ensure the algorithms run smoothly:

- Use matching IDs in both files.

- Set `observed_label` to 1 for selected 1-nodes and 0 for the others. The CSV reader preserves IDs as strings and labels as integers.

- Do not include self-loops or repeated edges.

### Trees or forests

As in the paper, we considered directed trees, i.e. edges point from parent to child. Replace the filenames below with your file paths:

``` r
source("load.R")
tree <- read_example_graph("tree_edges.csv", "tree_nodes.csv")
k <- sum(igraph::V(tree)$observed_label)  # Or choose another number of 1-nodes

optimal <- get_max_influence_and_label(tree, k, verbose = FALSE)
greedy <- maximize_influence_fixed_k_greedy(tree, k)
```

The scores are contained in `optimal$value` and `greedy$influence`, while their `label` vectors give the resulting 0/1 placements in node-file order.

### Networks

Use the same CSV columns for your directed network. In this case, nodes labeled 1 supply the initial seed set:

``` r
source("load.R")
network <- read_example_graph("network_edges.csv", "network_nodes.csv")
seeds <- which(igraph::V(network)$observed_label == 1)

set.seed(123)
simulation <- cascade_model(network, seeds, p = 0.3, nsim = 100, verbose = FALSE)
counts <- vapply(simulation$result, `[[`, integer(1), "n_infected")
mean(counts)
```

In this setting, `p` is the activation probability for each edge, and `nsim` is the number of repetitions.

## License and citation

The code and examples use the MIT License, which allows use,
modification, and redistribution, including commercial use, provided the copyright
and license notices are preserved. 

Please use the following citation when using these algorithms:

```
@misc{https://doi.org/10.48550/arxiv.2609.24398,
  doi = {10.48550/ARXIV.2609.24398},
  url = {https://arxiv.org/abs/2609.24398},
  author = {Di Marco,  Niccolò and Cinelli,  Matteo and Nakano,  Shinichi and Frosini,  Andrea},
  keywords = {Social and Information Networks (cs.SI),  Computers and Society (cs.CY),  FOS: Computer and information sciences},
  title = {Optimal and heuristic strategies for evaluating the influence of coordinated behavior in information cascades and retweet networks},
  publisher = {arXiv},
  year = {2026},
  copyright = {Creative Commons Attribution 4.0 International}
}
```
