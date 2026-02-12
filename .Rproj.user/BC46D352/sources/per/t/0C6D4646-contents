# ============================================================
# Script: analysis.R
# Purpose:
#   Reproduce analyses for inversion based phylogeny in cotton and radish
#   including NJ reconstruction, bipartition overlap, MAST, cophenetic correlation,
#   burden preserved randomization, bootstrap stability, and noise perturbation.
#
# Inputs:
#   data/tree3.csv   cotton inversion presence absence matrix
#   data/2tree.csv   radish inversion presence absence matrix
#
# Outputs:
#   figures/*.pdf
#   figures/*.csv
#   figures/*.txt
# ============================================================

library(ape)
library(phangorn)
library(ggplot2)
library(dplyr)

# ------------------------------
#functions
# ------------------------------

read_bin_matrix <- function(path) {
  mat <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  mat[is.na(mat)] <- 0
  mat[mat != 0] <- 1
  mat
}

nj_from_matrix <- function(mat) {
  phy <- phyDat(mat, type = "USER", levels = c(0, 1))
  NJ(dist.hamming(phy))
}

bipartitions_as_strings <- function(tree) {
  parts <- prop.part(tree)
  parts <- parts[sapply(parts, function(p) length(unlist(p)) < length(tree$tip.label))]
  unique(vapply(parts, function(p) {
    tips <- sort(tree$tip.label[unlist(p)])
    paste(tips, collapse = "|")
  }, character(1)))
}

safe_mast_size <- function(tree1, tree2) {
  tryCatch({
    m <- mast(tree1, tree2)
    length(m$tip.label)
  }, error = function(e) NA_integer_)
}

cophenetic_corr <- function(tree1, tree2) {
  common <- intersect(tree1$tip.label, tree2$tip.label)
  d1 <- cophenetic(keep.tip(tree1, common))
  d2 <- cophenetic(keep.tip(tree2, common))
  cor(as.vector(d1), as.vector(d2), method = "pearson", use = "complete.obs")
}

randomize_burden_preserved <- function(mat) {
  out <- matrix(0, nrow(mat), ncol(mat))
  for (i in seq_len(nrow(mat))) {
    ones <- sum(mat[i, ] == 1)
    if (ones > 0) out[i, sample.int(ncol(mat), ones)] <- 1
  }
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  out
}

add_noise_zeros_to_one <- function(mat, n) {
  zeros <- which(mat == 0)
  if (length(zeros) == 0) return(mat)
  n <- min(n, length(zeros))
  mat[sample(zeros, n)] <- 1
  mat
}

randomization_bip_overlap_thesis_style <- function(csv_path, n_iter = 100, seed = 123) {
  set.seed(seed)
  
  load_numeric <- function(filepath) {
    d <- read.csv(filepath, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    d <- as.data.frame(lapply(d, function(x) {
      x <- suppressWarnings(as.numeric(as.character(x)))
      x[is.na(x)] <- 0
      x
    }))
    d
  }
  
  clean <- function(mat) {
    mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
    mat <- mat[!duplicated(mat), , drop = FALSE]
    if (nrow(mat) < 4) return(NULL)
    repeat {
      phy <- try(phyDat(mat, type = "USER", levels = c(0, 1)), silent = TRUE)
      if (inherits(phy, "try-error")) return(NULL)
      dist_mat <- try(as.matrix(dist.hamming(phy)), silent = TRUE)
      if (inherits(dist_mat, "try-error")) return(NULL)
      na_rows <- unique(rownames(dist_mat)[rowSums(is.na(dist_mat)) > 0])
      if (length(na_rows) == 0) break
      mat <- mat[!rownames(mat) %in% na_rows, , drop = FALSE]
      if (nrow(mat) < 4) return(NULL)
    }
    mat
  }
  
  randomize_matrix <- function(data) {
    randomized <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    for (r in seq_len(nrow(data))) {
      ones <- sum(data[r, ] == 1)
      if (ones > 0 && ones <= ncol(data)) randomized[r, sample.int(ncol(data), ones)] <- 1
    }
    colnames(randomized) <- colnames(data)
    rownames(randomized) <- rownames(data)
    randomized
  }
  
  bipartitions_from_mat <- function(mat) {
    if (is.null(mat)) return(character(0))
    phy <- try(phyDat(mat, type = "USER", levels = c(0, 1)), silent = TRUE)
    if (inherits(phy, "try-error")) return(character(0))
    tree <- try(NJ(dist.hamming(phy)), silent = TRUE)
    if (inherits(tree, "try-error")) return(character(0))
    bip <- prop.part(tree)
    bip <- bip[sapply(bip, length) < length(tree$tip.label)]
    sapply(bip, function(part) paste(sort(tree$tip.label[unlist(part)]), collapse = "|"))
  }
  
  raw <- load_numeric(csv_path)
  orig_t <- clean(t(raw))
  bip_orig <- bipartitions_from_mat(orig_t)
  
  out <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    rand <- randomize_matrix(raw)
    rand_t <- clean(t(rand))
    bip_rand <- bipartitions_from_mat(rand_t)
    out[i] <- length(intersect(bip_orig, bip_rand))
  }
  
  out
}


# ------------------------------
# 1 Cotton NJ tree and reference tree
# ------------------------------

cot_mat <- read_bin_matrix("data/tree3.csv")
cot_tree <- nj_from_matrix(cot_mat)

cot_ref_newick <- "((Gm,Gt),((XLZ7,(ZMS24,NDM8)),(Gd,(Tanguis,(Xinhai21,(Hai7124,(Pima90,3-79)))))));"
cot_ref_tree <- read.tree(text = cot_ref_newick)

pdf("figures/cotton_nj_tree.pdf", width = 7, height = 6)
plot(cot_tree, main = "Cotton NJ tree from inversions")
dev.off()

pdf("figures/cotton_reference_tree.pdf", width = 7, height = 6)
plot(cot_ref_tree, main = "Cotton reference tree")
dev.off()

# ------------------------------
# 2 Cotton bipartition overlap
# ------------------------------

cot_bip <- bipartitions_as_strings(cot_tree)
cot_ref_bip <- bipartitions_as_strings(cot_ref_tree)
cot_shared <- intersect(cot_bip, cot_ref_bip)

writeLines(cot_shared, "figures/cotton_shared_bipartitions.txt")

cot_summary <- data.frame(
  dataset = "cotton",
  n_bip_nj = length(cot_bip),
  n_bip_ref = length(cot_ref_bip),
  n_shared = length(cot_shared)
)

print(cot_summary)

# ------------------------------
# 3 Cotton MAST and cophenetic correlation
# ------------------------------

cot_mast <- safe_mast_size(cot_tree, cot_ref_tree)
cot_cc <- cophenetic_corr(cot_tree, cot_ref_tree)

cat("Cotton MAST size:", cot_mast, "\n")
cat("Cotton cophenetic correlation:", cot_cc, "\n")

# ------------------------------
# 4 Burden preserved randomization test on bipartition overlap
# ------------------------------

set.seed(123)

n_iter <- 100
rand_shared_counts <- numeric(n_iter)

for (i in seq_len(n_iter)) {
  rand_mat <- randomize_burden_preserved(cot_mat)
  rand_tree <- nj_from_matrix(rand_mat)
  rand_bip <- bipartitions_as_strings(rand_tree)
  rand_shared_counts[i] <- length(intersect(cot_bip, rand_bip))
}

rand_df <- data.frame(shared_bipartitions = rand_shared_counts)

write.csv(rand_df, "figures/cotton_randomization_shared_bipartitions.csv", row.names = FALSE)

cat("Randomization mean:", mean(rand_shared_counts), "\n")
cat("Randomization sd:", sd(rand_shared_counts), "\n")
cat("Randomization min:", min(rand_shared_counts), "\n")
cat("Randomization max:", max(rand_shared_counts), "\n")

# ------------------------------
# 5 Bootstrap analysis
# ------------------------------

set.seed(42)

n_boot <- 100

boot_mast_vs_original <- numeric(n_boot)
boot_mast_vs_reference <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  sampled_cols <- sample.int(ncol(cot_mat), ncol(cot_mat), replace = TRUE)
  boot_mat <- cot_mat[, sampled_cols, drop = FALSE]
  boot_tree <- nj_from_matrix(boot_mat)
  
  boot_mast_vs_original[i] <- safe_mast_size(cot_tree, boot_tree)
  boot_mast_vs_reference[i] <- safe_mast_size(boot_tree, cot_ref_tree)
}

# Save results
boot_df <- data.frame(
  mast_vs_original = boot_mast_vs_original,
  mast_vs_reference = boot_mast_vs_reference
)

write.csv(boot_df, "figures/cotton_bootstrap_mast.csv", row.names = FALSE)

# Plot thesis metric: MAST vs reference
pdf("figures/cotton_bootstrap_mast_vs_reference_hist.pdf", width = 7, height = 5)
ggplot(data.frame(mast = boot_mast_vs_reference), aes(x = mast)) +
  geom_bar() +
  labs(title = "Cotton bootstrap MAST vs reference", x = "MAST size", y = "Count")
dev.off()

cat("Bootstrap (vs original) mean:", mean(boot_mast_vs_original), "\n")
cat("Bootstrap (vs reference) mean:", mean(boot_mast_vs_reference), "\n")
cat("Bootstrap (vs reference) min:", min(boot_mast_vs_reference), "\n")
cat("Bootstrap (vs reference) max:", max(boot_mast_vs_reference), "\n")

# ------------------------------
# 6 Noise perturbation by flipping zeros to ones
# ------------------------------

set.seed(42)

noise_levels <- seq(0, sum(cot_mat == 0), by = 30)
n_reps <- 10

mast_data <- data.frame()

for (noise_n in noise_levels) {
  for (rep in seq_len(n_reps)) {
    noisy_mat <- add_noise_zeros_to_one(cot_mat, noise_n)
    noisy_tree <- nj_from_matrix(noisy_mat)
    mast_val <- safe_mast_size(cot_tree, noisy_tree)
    mast_data <- rbind(mast_data, data.frame(noise = noise_n, rep = rep, mast = mast_val))
  }
  cat("Finished noise level:", noise_n, "\n")
}

mast_summary <- mast_data %>%
  group_by(noise) %>%
  summarise(mean_mast = mean(mast, na.rm = TRUE), sd = sd(mast, na.rm = TRUE), .groups = "drop")

write.csv(mast_summary, "figures/cotton_noise_mast_summary.csv", row.names = FALSE)

pdf("figures/cotton_noise_mast_curve.pdf", width = 7, height = 5)
ggplot(mast_summary, aes(x = noise, y = mean_mast)) +
  geom_line() +
  geom_point() +
  labs(title = "Cotton noise perturbation", x = "Number of zeros flipped to one", y = "Mean MAST size")
dev.off()



# ============================================================
# Radish dataset
# ============================================================

# ------------------------------
# Radish 1 NJ reconstruction and reference tree
# ------------------------------

rad_mat <- read_bin_matrix("data/2tree.csv")
rad_tree <- nj_from_matrix(rad_mat)

# Reference tree from Zhang et al. 2021 (as used in rs001.Rmd)
rad_ref_newick <- "((((RS00,RS04),RS02),((RS05,RS07),RS01)),(((RS06,RS10),RS03),(RS08,RS09)));"
rad_ref_tree <- read.tree(text = rad_ref_newick)

pdf("figures/radish_nj_tree.pdf", width = 7, height = 6)
plot(rad_tree, main = "Radish NJ tree from inversions")
dev.off()

pdf("figures/radish_reference_tree.pdf", width = 7, height = 6)
plot(rad_ref_tree, main = "Radish reference tree (Zhang et al. 2021)")
dev.off()

# ------------------------------
# Radish 2 Bipartition overlap
# ------------------------------

rad_bip <- bipartitions_as_strings(rad_tree)
rad_ref_bip <- bipartitions_as_strings(rad_ref_tree)
rad_shared <- intersect(rad_bip, rad_ref_bip)

writeLines(rad_shared, "figures/radish_shared_bipartitions.txt")

rad_summary <- data.frame(
  dataset = "radish",
  n_bip_nj = length(rad_bip),
  n_bip_ref = length(rad_ref_bip),
  n_shared = length(rad_shared)
)
print(rad_summary)

# ------------------------------
# Radish 3 MAST and cophenetic correlation
# ------------------------------

rad_mast <- safe_mast_size(rad_tree, rad_ref_tree)
rad_cc <- cophenetic_corr(rad_tree, rad_ref_tree)

cat("Radish MAST size:", rad_mast, "\n")
cat("Radish cophenetic correlation:", rad_cc, "\n")

# ------------------------------
# Radish 4 Bipartition randomization test (thesis style)
# Row-burden preserved shuffling + transpose/clean pipeline used for the thesis table
# ------------------------------

rad_rand_shared <- randomization_bip_overlap_thesis_style(
  csv_path = "data/2tree.csv",
  n_iter = 100,
  seed = 123
)

write.csv(
  data.frame(shared_bipartitions = rad_rand_shared),
  "figures/radish_randomization_shared_bipartitions.csv",
  row.names = FALSE
)

cat("Radish randomization mean:", mean(rad_rand_shared), "\n")
cat("Radish randomization sd:", sd(rad_rand_shared), "\n")
cat("Radish randomization min:", min(rad_rand_shared), "\n")
cat("Radish randomization max:", max(rad_rand_shared), "\n")

# ------------------------------
# Radish 5 Bootstrap (vs reference)
# Thesis style comparison: bootstrap replicate trees vs published reference
# ------------------------------

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(42)

n_boot <- 100
rad_boot_mast_ref <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  sampled_cols <- sample.int(ncol(rad_mat), ncol(rad_mat), replace = TRUE)
  boot_mat <- rad_mat[, sampled_cols, drop = FALSE]
  boot_tree <- nj_from_matrix(boot_mat)
  rad_boot_mast_ref[i] <- safe_mast_size(boot_tree, rad_ref_tree)
}

write.csv(
  data.frame(mast_vs_reference = rad_boot_mast_ref),
  "figures/radish_bootstrap_mast_vs_reference.csv",
  row.names = FALSE
)

cat("Radish bootstrap mean vs reference:", mean(rad_boot_mast_ref), "\n")
cat("Radish bootstrap min vs reference:", min(rad_boot_mast_ref), "\n")
cat("Radish bootstrap max vs reference:", max(rad_boot_mast_ref), "\n")

pdf("figures/radish_bootstrap_mast_vs_reference_hist.pdf", width = 7, height = 5)
ggplot(data.frame(mast = rad_boot_mast_ref), aes(x = mast)) +
  geom_bar() +
  labs(title = "Radish bootstrap MAST vs reference", x = "MAST size", y = "Count")
dev.off()
