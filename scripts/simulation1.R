# ============================================================
# Simulation 1
#   Evaluate how binary tip states (Dollo-like single gain + inheritance)
#   respond to increasing perturbation (tip flipping), measured by:
#     1) MAST size between original vs perturbed NJ trees
#     2) Cophenetic correlation between original vs perturbed NJ trees
#     3) Local stability proxy: preserved sister pairs (cherries)
#
# Outputs:
#   outputs/sim1_tree_log.csv
#   outputs/sim1_rep_results.csv
#   outputs/sim1_tree_summary.csv
#   outputs/sim1_bin_summary.csv
#   outputs/sim1_correlation_mast_cherries.csv
#   figures/sim1_mast_all_bins.pdf
#   figures/sim1_cophenetic_all_bins.pdf
#   figures/sim1_correlation_all_bins.pdf
# ============================================================

library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(42)

dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ------------------------------
# Key parameters (match thesis)
# ------------------------------
ntips <- 32
n_bins <- 20
trees_per_bin <- 20
n_reps <- 15
max_attempts <- 20000

# ------------------------------
# Bin definitions for initial prop(1)
# 20 bins of width 5%
# ------------------------------
bin_edges <- seq(0, 1, length.out = n_bins + 1)
bin_labels <- paste0(
  sprintf("%02.0f", bin_edges[-length(bin_edges)] * 100),
  "-",
  sprintf("%02.0f", bin_edges[-1] * 100),
  "%"
)

get_bin_index <- function(p) {
  idx <- which(p >= bin_edges[-length(bin_edges)] & p < bin_edges[-1])
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

# ------------------------------
# Dollo-like propagation:
# one internal node gains state 1, and all its descendants inherit it
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

# ------------------------------
# Reconstruct NJ tree from 1 binary character (tip states)
# Using Manhattan distance on a 1-column matrix (equivalent to Hamming here)
# ------------------------------
nj_from_tip_states <- function(tip_states, tip_labels) {
  char_mat <- matrix(tip_states, ncol = 1)
  rownames(char_mat) <- tip_labels
  nj(dist(char_mat, method = "manhattan"))
}

mast_size <- function(t1, t2) {
  Ntip(mast(unroot(t1), unroot(t2)))
}

cophenetic_corr <- function(t1, t2) {
  common <- intersect(t1$tip.label, t2$tip.label)
  d1 <- cophenetic(t1)[common, common]
  d2 <- cophenetic(t2)[common, common]
  cor(as.vector(d1), as.vector(d2), use = "complete.obs")
}

# ------------------------------
# Sister pairs (cherries) as local topology units
# Return as canonical strings: "tipA|tipB"
# ------------------------------
get_cherries <- function(tr) {
  ed <- tr$edge
  kids <- split(ed[, 2], ed[, 1])
  cherries <- lapply(kids, function(v) {
    tips <- v[v <= Ntip(tr)]
    if (length(tips) == 2) tr$tip.label[tips] else NULL
  })
  cherries <- Filter(Negate(is.null), cherries)
  unique(vapply(cherries, function(x) paste(sort(x), collapse = "|"), character(1)))
}

# ------------------------------
# Storage
# ------------------------------
bin_counts <- rep(0, n_bins)
names(bin_counts) <- bin_labels

tree_log <- data.frame(
  bin = character(),
  tree_in_bin = integer(),
  prop1 = numeric(),
  stringsAsFactors = FALSE
)

rep_results <- data.frame(
  bin = character(),
  tree_id = integer(),
  tree_in_bin = integer(),
  prop1 = numeric(),
  k = integer(),
  rep = integer(),
  mast = integer(),
  coph = numeric(),
  preserved_cherries = integer(),
  stringsAsFactors = FALSE
)

# ============================================================
# Main loop:
# collect 20 trees per bin (if possible), then perturb each tree
# ============================================================
attempts <- 0

while (any(bin_counts < trees_per_bin) && attempts < max_attempts) {
  attempts <- attempts + 1
  
  # 1) Random unrooted tree
  tree <- rtree(ntips)
  
  # 2) Assign Dollo-like binary states on all nodes
  tot <- ntips + tree$Nnode
  children <- split(tree$edge[, 2], tree$edge[, 1])
  
  gain_node <- sample((ntips + 1):tot, 1)
  states <- integer(tot)
  states[gain_node] <- 1
  states <- propagate_down(children, states, gain_node)
  
  tip_states <- states[1:ntips]
  prop1 <- mean(tip_states)
  
  # 3) Bin by initial prop(1)
  bin_id <- get_bin_index(prop1)
  if (is.na(bin_id)) next
  if (bin_counts[bin_id] >= trees_per_bin) next
  
  # Register tree
  bin_counts[bin_id] <- bin_counts[bin_id] + 1
  bin_label <- bin_labels[bin_id]
  tree_in_bin <- bin_counts[bin_id]
  tree_id <- sum(bin_counts)
  
  tree_log <- rbind(tree_log, data.frame(
    bin = bin_label,
    tree_in_bin = tree_in_bin,
    prop1 = prop1
  ))
  
  cat(sprintf("Collected tree %2d in bin %s (prop1 = %.3f)\n",
              tree_in_bin, bin_label, prop1))
  
  # 4) Reference NJ tree from original tip states
  ref_tree <- nj_from_tip_states(tip_states, tree$tip.label)
  ref_cherries <- get_cherries(unroot(ref_tree))
  
  # 5) Progressive perturbation: flip k tips, repeat 15 times
  for (k in 1:ntips) {
    for (rr in 1:n_reps) {
      flip_idx <- sample(seq_len(ntips), k, replace = FALSE)
      tip_states2 <- tip_states
      tip_states2[flip_idx] <- 1 - tip_states2[flip_idx]
      
      flip_tree <- nj_from_tip_states(tip_states2, tree$tip.label)
      
      mast_val <- mast_size(ref_tree, flip_tree)
      coph_val <- cophenetic_corr(ref_tree, flip_tree)
      
      flip_cherries <- get_cherries(unroot(flip_tree))
      preserved <- sum(ref_cherries %in% flip_cherries)
      
      rep_results <- rbind(rep_results, data.frame(
        bin = bin_label,
        tree_id = tree_id,
        tree_in_bin = tree_in_bin,
        prop1 = prop1,
        k = k,
        rep = rr,
        mast = mast_val,
        coph = coph_val,
        preserved_cherries = preserved
      ))
    }
  }
}

# ------------------------------
# Save raw results
# ------------------------------
write.csv(tree_log, "outputs/sim1_tree_log.csv", row.names = FALSE)
write.csv(rep_results, "outputs/sim1_rep_results.csv", row.names = FALSE)
cat("Saved: outputs/sim1_tree_log.csv\n")
cat("Saved: outputs/sim1_rep_results.csv\n")

# ------------------------------
# Summaries:
# tree_summary: each tree gives one curve (mean over 15 reps at each k)
# bin_summary: mean curve per bin (average across 20 trees)
# ------------------------------
tree_summary <- rep_results %>%
  group_by(bin, tree_id, tree_in_bin, prop1, k) %>%
  summarise(
    mean_mast = mean(mast, na.rm = TRUE),
    mean_coph = mean(coph, na.rm = TRUE),
    mean_preserved_cherries = mean(preserved_cherries, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(tree_summary, "outputs/sim1_tree_summary.csv", row.names = FALSE)
cat("Saved: outputs/sim1_tree_summary.csv\n")

bin_summary <- tree_summary %>%
  group_by(bin, k) %>%
  summarise(
    mean_mast = mean(mean_mast, na.rm = TRUE),
    sd_mast = sd(mean_mast, na.rm = TRUE),
    mean_coph = mean(mean_coph, na.rm = TRUE),
    sd_coph = sd(mean_coph, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(bin_summary, "outputs/sim1_bin_summary.csv", row.names = FALSE)
cat("Saved: outputs/sim1_bin_summary.csv\n")

# ------------------------------
# Plot: 20 lines per bin panel (match thesis figure)
# ------------------------------
tree_summary$bin <- factor(tree_summary$bin, levels = bin_labels)

p_mast_bins_20lines <- ggplot(
  tree_summary,
  aes(x = k, y = mean_mast, group = interaction(bin, tree_id))
) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  facet_wrap(~ bin, nrow = 5) +
  labs(
    title = "MAST sensitivity across 20 bins of prop1",
    x = "Number of Flipped Tips",
    y = "Mean MAST"
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/sim1_mast_all_bins.pdf", p_mast_bins_20lines, width = 12, height = 10)
cat("Saved: figures/sim1_mast_all_bins.pdf\n")

p_coph_bins_20lines <- ggplot(
  tree_summary,
  aes(x = k, y = mean_coph, group = interaction(bin, tree_id))
) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  facet_wrap(~ bin, nrow = 5) +
  labs(
    title = "Cophenetic correlation across 20 bins of prop1",
    x = "Number of Flipped Tips",
    y = "Mean cophenetic correlation"
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/sim1_cophenetic_all_bins.pdf", p_coph_bins_20lines, width = 12, height = 10)
cat("Saved: figures/sim1_cophenetic_all_bins.pdf\n")

# ------------------------------
# Correlation: preserved cherries vs mean MAST
# Compute correlation across trees within each bin at each k
# Some k can be undefined (sd = 0), reported as NA
# ------------------------------
cor_results <- tree_summary %>%
  group_by(bin, k) %>%
  summarise(
    cor_mast_cherries = {
      x <- mean_mast
      y <- mean_preserved_cherries
      if (length(x) < 3 || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
        NA_real_
      } else {
        suppressWarnings(cor(x, y, use = "complete.obs"))
      }
    },
    .groups = "drop"
  )

write.csv(cor_results, "outputs/sim1_correlation_mast_cherries.csv", row.names = FALSE)
cat("Saved: outputs/sim1_correlation_mast_cherries.csv\n")

cor_results$bin <- factor(cor_results$bin, levels = bin_labels)
cor_plot_data <- cor_results %>% filter(!is.na(cor_mast_cherries))

plot_list <- lapply(levels(cor_plot_data$bin), function(bb) {
  dfb <- cor_plot_data %>% filter(bin == bb)
  ggplot(dfb, aes(x = k, y = cor_mast_cherries)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 0.8) +
    labs(title = paste("Bin", bb), x = "Flipped tips (k)", y = "Correlation") +
    ylim(-1, 1) +
    theme_minimal(base_size = 10)
})

big_plot <- wrap_plots(plotlist = plot_list, ncol = 4)
ggsave("figures/sim1_correlation_all_bins.pdf", big_plot, width = 16, height = 18)
cat("Saved: figures/sim1_correlation_all_bins.pdf\n")
