library(DEP)
library(tidyverse)
library(SummarizedExperiment)
library(openxlsx)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(stringr)

# Parameters
max_pep_thresh <- 3
logFC_iact     <- 1
p_val_thresh   <- 0.01
max_NA         <- 7

# Load data
pulldown_combi <- read.xlsx(file.path("2021-143_MaxLFQ_Intensities_Th17_TReg_Pulldown_together_MBR_nocrossnorm_Core4_UD.xlsx"), sheet = 3, colNames = FALSE)
names(pulldown_combi) <- make.names(pulldown_combi[2,], unique = TRUE)
pulldown_combi <- pulldown_combi %>% filter(Protein.Group != "Protein.Group")
pulldown_combi$prot.names <- str_remove(pulldown_combi$prot.names, "_MOUSE")

pulldown_combi_base <- read.xlsx(file.path(data.folder, "2021-143_MaxLFQ_Intensities_Th17_TReg_Pulldown_together_MBR_nocrossnorm_Core4_UD.xlsx"), sheet = 2) %>%
  select(-X1) %>%
  mutate(prot.names = str_remove(prot.names, "_MOUSE"))

experiments_combi <- grep("Exp", colnames(pulldown_combi), value = TRUE)
anno_combi <- data.frame(
  row.names = experiments_combi,
  Type = str_extract(experiments_combi, "WT|Bio"),
  Exp  = str_extract(experiments_combi, "Exp..[0-9]*")
)
anno_combi$Cell <- factor(rep(c("Treg", "Th17"), each = length(experiments_combi) / 2), levels = c("Th17", "Treg"))

experiments_Th17 <- rownames(anno_combi %>% filter(Cell == "Th17"))
pulldown_Th17 <- pulldown_combi %>%
  select(prot.names, Protein.Group, all_of(experiments_Th17), `SIGNIFICANT_BH...Log2.Filter.1`) %>%
  filter(`SIGNIFICANT_BH...Log2.Filter.1` != "ns")

pulldown_Th17_base <- pulldown_combi_base %>%
  select(prot.names, Protein.Group, max_pep, all_of(experiments_Th17))

anno_Th17 <- data.frame(
  row.names = experiments_Th17,
  Type = str_extract(experiments_Th17, "WT|Bio"),
  Exp  = str_extract(experiments_Th17, "Exp..[0-9]*")
)

# Design and contrast
design_th17 <- model.matrix(~0 + Type + Exp, data = anno_Th17)
contr.matrix <- makeContrasts(BiovsWT = TypeBio - TypeWT, levels = design_th17)

# Input matrix
limma_th17_input_raw <- log2(pulldown_Th17_base[experiments_Th17])
rownames(limma_th17_input_raw) <- pulldown_Th17_base$prot.names
limma_th17_input_raw <- limma_th17_input_raw[pulldown_Th17_base$max_pep >= max_pep_thresh, ]

# DEP formatting
experimental_design_th17 <- data.frame(
  label = experiments_Th17,
  condition = str_extract(experiments_Th17, "WT|Bio"),
  technical_replicate = str_extract(experiments_Th17, "[1-3]$"),
  Exp = str_extract(experiments_Th17, "Exp..[0-9][0-9]")
) %>%
  mutate(replicate = paste(Exp, technical_replicate, sep = "_"))

bool_col <- which(colnames(pulldown_Th17_base) %in% experiments_Th17)
unique_prot <- make_unique(pulldown_Th17_base, "prot.names", "Protein.Group", delim = ";")
dep_input_raw <- make_se(unique_prot[unique_prot$max_pep >= max_pep_thresh, ], bool_col, experimental_design_th17)

# MNAR detection and mixed imputation
proteins_MNAR <- get_df_long(dep_input_raw) %>%
  group_by(name, condition, Exp) %>%
  summarize(NAs = all(is.na(intensity)), .groups = "drop") %>%
  filter(NAs) %>%
  pull(name) %>%
  unique()

MNAR <- names(dep_input_raw) %in% proteins_MNAR

dep_input_raw_imp <- impute(
  dep_input_raw,
  fun = "mixed",
  randna = !MNAR,
  mar = "knn",
  mnar = "min"
)

limma_th17_input <- assay(dep_input_raw_imp)
colnames(limma_th17_input) <- colnames(limma_th17_input_raw)
rownames(limma_th17_input) <- rownames(limma_th17_input_raw)

# Limma fit
fit <- lmFit(limma_th17_input, design_th17)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

results_paired_th17 <- as.data.frame(topTable(fit, number = nrow(limma_th17_input)))
results_paired_th17_comp <- results_paired_th17

# Filtering based on max_NA threshold
bio_cols <- anno_Th17$Type == "Bio"
wt_cols  <- anno_Th17$Type == "WT"

invalid_bio <- rowSums(is.na(limma_th17_input_raw[, bio_cols])) > max_NA
invalid_wt  <- rowSums(is.na(limma_th17_input_raw[, wt_cols]))  > max_NA
valid_proteins <- !(invalid_bio | invalid_wt)

results <- results_paired_th17[rownames(limma_th17_input_raw)[valid_proteins], ]
results <- results %>% filter(adj.P.Val < p_val_thresh, logFC >= logFC_iact)

write.csv(results,"results_limma_th17_vs_wt.csv")
