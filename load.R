# Load the algorithm functions with source("load.R").
.code_stack <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)
}, character(1))
.code_file <- tail(.code_stack[nzchar(.code_stack)], 1)
if (!length(.code_file)) stop("Load this file with source().")
.code_dir <- dirname(normalizePath(.code_file))
for (.code_name in c("influence_constrained_k_binary_tree.R",
                     "get_max_influence_and_label_forest.R",
                     "get_max_influence_and_label.R",
                     "influence_constrained_k_greedy.R", "cascade_model.R",
                     "utils.R", "plots.R")) {
  sys.source(file.path(.code_dir, "src", .code_name), envir = environment())
}
rm(.code_stack, .code_file, .code_dir, .code_name)
