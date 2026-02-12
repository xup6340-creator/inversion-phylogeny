# ============================================================
# Simulation 2 
# Controlled initial trait-1 proportion (10–50%), then tip flipping
# Output:
#   figures/sim2_example_trees.pdf   (5 example trees with red/blue tips)
#   figures/sim2_mast_curve.pdf      (5-group mean MAST curve, like sim2fin.png)
#   outputs/sim2_tree_log.csv
#   outputs/sim2_results.csv
#   outputs/sim2_summary.csv
# ============================================================

library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(ggtree)

# ------------------------------
# Reproducibility and folders
# ------------------------------
set.seed(42)
dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ------------------------------
# Parameters
# ------------------------------
ntips <- 32
n_reps <- 15

# 5 bins for controlled prop1 (as in thesis)
bins <- list(
  c(0.0, 0.1),
  c(0.1, 0.2),
  c(0.2, 0.3),
  c(0.3, 0.4),
  c(0.4, 0.5)
)
bin_labels <- c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%")

# How close we require the achieved proportion to be
# In practice, exact matching is impossible because prop1 jumps by 1/32.
# This tolerance keeps selection stable.
tol <- 0.02

# Safety limit to avoid infinite loops
max_attempts_per_bin <- 20000

# ------------------------------
# Helper functions
# ------------------------------
propagate_down <- function(children, states, node) {
  kids <- children[[as.character(node)]]
  if (!is.null(kids)) {
    for (ch in kids) {
      states[ch] <- 1
      states <- propagate_down(children, states, ch)
    }
  }
  states
}

nj_from_tip_states <- function(tip_states, tip_labels) {
  mat <- matrix(tip_states, ncol = 1)
  rownames(mat) <- tip_labels
  nj(dist(mat, method = "manhattan"))
}

mast_size <- function(t1, t2) {
  Ntip(mast(unroot(t1), unroot(t2)))
}

# Plot one tree with colored tip points (red=1, blue=0)
plot_tree_with_states <- function(tree, tip_states, main = "") {
  plot(tree, show.tip.label = FALSE, main = main)
  tipcols <- ifelse(tip_states == 1, "red", "blue")
  tiplabels(pch = 19, col = tipcols, cex = 1)
}

# ------------------------------
# 1) Pick one representative tree per bin (controlled prop1)
# ------------------------------
tree_log <- data.frame(
  bin = character(),
  prop1 = numeric(),
  ones = integer(),
  zeros = integer(),
  gain_node = integer(),
  attempts = integer(),
  stringsAsFactors = FALSE
)

picked <- vector("list", length(bins))

for (b in seq_along(bins)) {
  target_mid <- mean(bins[[b]])
  attempts <- 0
  
  repeat {
    attempts <- attempts + 1
    if (attempts > max_attempts_per_bin) {
      stop(paste0(
        "Failed to find a tree for bin ", bin_labels[b],
        " within max attempts. Increase tol or max_attempts_per_bin."
      ))
    }
    
    tree <- rtree(ntips)
    tot <- ntips + tree$Nnode
    children <- split(tree$edge[, 2], tree$edge[, 1])
    
    gain_node <- sample((ntips + 1):tot, 1)
    states <- integer(tot)
    states[gain_node] <- 1
    states <- propagate_down(children, states, gain_node)
    
    tip_states <- states[1:ntips]
    prop1 <- mean(tip_states)
    
    if (abs(prop1 - target_mid) <= tol &&
        prop1 > bins[[b]][1] && prop1 <= bins[[b]][2]) {
      
      picked[[b]] <- list(
        tree = tree,
        tip_states = tip_states,
        prop1 = prop1,
        gain_node = gain_node
      )
      
      ones <- sum(tip_states)
      zeros <- ntips - ones
      
      tree_log <- rbind(tree_log, data.frame(
        bin = bin_labels[b],
        prop1 = prop1,
        ones = ones,
        zeros = zeros,
        gain_node = gain_node,
        attempts = attempts
      ))
      
      cat(sprintf(
        "Selected %s: prop1=%.3f (ones=%d, zeros=%d) after %d attempts\n",
        bin_labels[b], prop1, ones, zeros, attempts
      ))
      
      break
    }
  }
}

write.csv(tree_log, "outputs/sim2_tree_log.csv", row.names = FALSE)
cat("Saved: outputs/sim2_tree_log.csv\n")

# ------------------------------
# 2) Save example trees with node + tip states (ggtree style)
#    This matches your screenshot: every node has a red/blue marker.
# ------------------------------

plot_tree_with_node_states <- function(tree, states, ntips, title_txt = "") {
  tot <- ntips + tree$Nnode
  df_states <- data.frame(
    node = 1:tot,
    state = factor(states, levels = c(0, 1))
  )
  
  ggtree(tree, layout = "rectangular") %<+% df_states +
    geom_nodepoint(aes(color = state), size = 2.8) +
    geom_tippoint(aes(color = state), size = 2.8) +
    scale_color_manual(values = c("0" = "blue", "1" = "red")) +
    ggtitle(title_txt) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

pdf("figures/sim2_example_trees.pdf", width = 12, height = 8)
for (b in seq_along(bins)) {
  obj <- picked[[b]]
  
  # Reconstruct full node states for plotting (tips + internal nodes)
  # We stored only tip_states in picked; we must regenerate node states.
  # Easiest: re-run gain on the same tree using the stored gain_node.
  tree <- obj$tree
  tot <- ntips + tree$Nnode
  children <- split(tree$edge[, 2], tree$edge[, 1])
  
  states <- integer(tot)
  states[obj$gain_node] <- 1
  states <- propagate_down(children, states, obj$gain_node)
  
  main_txt <- sprintf(
    "%d: %.1f%% red tips (%d tips), gain node = %d",
    b, mean(states[1:ntips]) * 100, sum(states[1:ntips]), obj$gain_node
  )
  
  p_state <- plot_tree_with_node_states(tree, states, ntips, main_txt)
  print(p_state)
}
dev.off()

cat("Saved: figures/sim2_example_trees.pdf\n")


# ------------------------------
# 3) Flipping experiment per bin (like Simulation 1, but only 5 groups)
# ------------------------------
results <- data.frame(
  bin = character(),
  k = integer(),
  rep = integer(),
  mast = integer(),
  stringsAsFactors = FALSE
)

for (b in seq_along(bins)) {
  obj <- picked[[b]]
  ref_tree <- nj_from_tip_states(obj$tip_states, obj$tree$tip.label)
  
  for (k in 1:ntips) {
    for (rr in 1:n_reps) {
      tip2 <- obj$tip_states
      flip_idx <- sample(seq_len(ntips), k, replace = FALSE)
      tip2[flip_idx] <- 1 - tip2[flip_idx]
      
      flip_tree <- nj_from_tip_states(tip2, obj$tree$tip.label)
      m <- mast_size(ref_tree, flip_tree)
      
      results <- rbind(results, data.frame(
        bin = bin_labels[b],
        k = k,
        rep = rr,
        mast = m
      ))
    }
  }
  
  cat(sprintf("Finished flipping: %s\n", bin_labels[b]))
}

write.csv(results, "outputs/sim2_results.csv", row.names = FALSE)
cat("Saved: outputs/sim2_results.csv\n")

# ------------------------------
# 4) Summarize and plot (one line per bin)
# ------------------------------
summary_df <- results %>%
  group_by(bin, k) %>%
  summarise(
    mean_mast = mean(mast, na.rm = TRUE),
    sd_mast = sd(mast, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, "outputs/sim2_summary.csv", row.names = FALSE)
cat("Saved: outputs/sim2_summary.csv\n")

summary_df$bin <- factor(summary_df$bin, levels = bin_labels)

p <- ggplot(summary_df, aes(x = k, y = mean_mast, color = bin)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  labs(
    title = "Simulation 2: Mean MAST under increasing tip flipping",
    x = "Number of Flipped Tips",
    y = "Mean MAST",
    color = "Trait 1 Ratio"
  ) +
  theme_minimal(base_size = 14)

ggsave("figures/sim2_mast_curve.pdf", p, width = 9, height = 5.5)
cat("Saved: figures/sim2_mast_curve.pdf\n")

print(p)
