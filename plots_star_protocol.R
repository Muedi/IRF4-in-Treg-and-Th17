# ======================================================
# === Install packages if not already installed
# ======================================================

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
  install.packages("ggseqlogo")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# ======================================================
# === Load required packages
# ======================================================

library(Biostrings)
library(data.table)
library(ggseqlogo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(readr)
library(ggrepel)

# ======================================================
# === 1) Heatmap: Double Motif Co-Occurrence
# ======================================================

setwd("~/phd/star/v3")

motif_pairs <- read.csv("combined_motifs_stats_5bp_csv.csv")

df <- motif_pairs %>%
  filter(Pulldown. == "both") %>%
  mutate(Norm_Count = Count_Th17 / max(Count_Th17, na.rm = TRUE)) %>%
  select(TF1, TF2, Norm_Count)

df_wide <- df %>%
  pivot_wider(names_from = TF2, values_from = Norm_Count, values_fill = 0)

row_sums <- rowSums(df_wide[,-1])
col_sums <- colSums(df_wide[,-1])

row_order <- order(-row_sums)
col_order <- names(sort(-col_sums))

df_sorted <- df_wide[row_order, c("TF1", col_order)]

mat <- as.matrix(df_sorted[,-1])
rownames(mat) <- df_sorted$TF1

pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
  main = "Double Motif Co-Occurrence (Th17, Pulldown: both)",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45
)

# ======================================================
# === 2) Volcano Plot: LIMMA Results
# ======================================================

results <- read_csv("results_limma_th17_vs_wt.csv")

results$Significance <- "Not Significant"
results$Significance[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Up"
results$Significance[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Down"

TFs <- unique(df_sorted$TF1)
TFs <- c(TFs, "EP300")

results$Label <- ifelse(
  results$...1 %in% TFs & results$Significance != "Not Significant",
  results$...1,
  NA
)

p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), size = 1.5) +
  geom_text_repel(
    aes(label = Label),
    size = 3,
    na.rm = TRUE,
    color = "black",
    segment.color = "black",
    segment.size = 0.3
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "Not Significant" = "grey")
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot with Significant TF1 Labels",
    x = "log2 Fold Change (Th1 vs. WT)",
    y = "-log10 Adjusted P-Value",
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

print(p)

# ======================================================
# === 3) Extract PWM from MATCH OUTPUT & PLOT LOGOS
# ======================================================

# Load your motif matches
match_data <- read.table("match_outpuz.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
colnames(match_data) <- c("ID", "Name", "Value1", "Strand", "Value2", "Value3", "Sequence", "Gene")

# Make sequences uppercase
match_data$Sequence <- toupper(match_data$Sequence)

# Get unique TFs
unique_motifs <- unique(match_data$Name)

# Loop over all motifs
for (motif in unique_motifs) {
  # Extract sequences for this motif
  seqs <- match_data[match_data$Name == motif, "Sequence"]
  seq_set <- DNAStringSet(seqs)
  
  # Make consensus matrix
  pwm <- consensusMatrix(seq_set, as.prob = FALSE, baseOnly = TRUE)
  pwm <- pwm[c("A", "C", "G", "T"), ]
  colnames(pwm) <- 1:ncol(pwm)
  
  # Save PWM file
  pwm_file <- paste0(motif, ".pwm")
  write.table(pwm, file = pwm_file, sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
  
  # Plot logo from PWM
  pwm_matrix <- as.matrix(pwm)
  logo <- ggseqlogo(pwm_matrix, method = "bits") +
    ggtitle(motif) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank())
  
  ggsave(paste0(motif, "_logo.png"), logo, width = 6, height = 2)
}

