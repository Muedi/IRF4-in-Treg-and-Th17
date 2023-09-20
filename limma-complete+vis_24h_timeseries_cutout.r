######################################################################################################
# Script: Limma Analysis of 24h Proteome and IRF4-Pulldown (Interactome) Data
# Author: Maximilian Sprang, Muedi
# Date: 20.09.2023
# Description: 
# This script Uses Limma to obtain information on differentially expressed proteins, from MS-Proteomic data. 
# In combination with a pull down, we can observe an Interactome, when we compare to non-pulled-down proteomes.
# DEPs imputation was ustilized.
# Multiple Visualizations are included below.
######################################################################################################
# libs
library(DEP)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
library(ggVennDiagram)
library(openxlsx)
library(limma)
library(EnhancedVolcano)
# library(diann)
library(openxlsx)
library(pheatmap)
library(stringr)
library(FactoMineR)
library(ggrepel)
library(viridis)

outdir <- "output_new/mouse_24h" #_max_pep-3_run_filter_6"
dir.create(outdir)

max_pep_thresh_prot = 2
max_pep_thresh = 3
max_NAs = 3

p_val_thresh = 0.01
logFC_prot = 0.5
logFC_iact = 1

################################################################################
################################################################################
#                                  Proteome                                    #      
################################################################################
################################################################################
# missing data mixed imp for Proteome limma results with DEP imputation.
# Proteome is filtered for:
# maxpep >= 2 
# runfilter: 
# th17_runfilter <- c(limma_wholeProt_input_raw[anno_wholeProt$Cell == "Th17"] %>% apply(1, count_na)) > 1
# treg_runfilter <- c(limma_wholeProt_input_raw[anno_wholeProt$Cell == "Treg"] %>% apply(1, count_na)) > 1
# results_prot <- results_prot[rownames(limma_wholeProt_input_raw[!(th17_runfilter & treg_runfilter),]), ] # remove Proteins, that do not have at least one condition with more than 8 non-NA samples 
# this means that we generate a boolean vector for both Th17 and Tregs
# holding a TRUE evertime a Protein has more than 1 NA in either of the groups.
# The two vectors are combined by a logical & and negated resulting in a dataset
# with only Proteins holdiung more than 8 non-NA values in the raw intensities.
# adj.P <= 0.05
# logFC > 0.5

folder.prot <- file.path("data", "time_series", "Proteom")
file.prot <- "2022-144_Proteom_Idefix_MaxLFQ_Intensities_AG.xlsx"

proteome_Th17vTreg_base <- read.xlsx(file.path(folder.prot, file.prot), sheet = 4)
proteome_Th17vTreg_base <- proteome_Th17vTreg_base %>% dplyr::select(-X1)
proteome_Th17vTreg_base$prot.names <- str_remove_all(proteome_Th17vTreg_base$prot.names, "_MOUSE")
proteome_Th17vTreg_base$max_pep <- as.numeric(proteome_Th17vTreg_base$max_pep)
# proteome_Th17vTreg_base[is.na(proteome_Th17vTreg_base)] <- 0
experiments_wholeProt <- colnames(proteome_Th17vTreg_base)
experiments_wholeProt <- experiments_wholeProt[grepl("Exp", experiments_wholeProt)]
experiments_wholeProt <- experiments_wholeProt[grepl("24h", experiments_wholeProt)] # only 24h

anno_wholeProt <- data.frame(row.names=experiments_wholeProt,"Type"=str_extract(experiments_wholeProt, "WT|Bio"))
celltypes_prot <- str_extract(experiments_wholeProt, "Th17|Treg")
anno_wholeProt["Cell"] <- factor(celltypes_prot, levels = c("Th17", "Treg"))
anno_wholeProt_row <- data.frame(row.names =proteome_Th17vTreg_base$prot.names,"max_pep"=proteome_Th17vTreg_base$max_pep)

# make input_raw
limma_wholeProt_input_raw <- log2(proteome_Th17vTreg_base[experiments_wholeProt])
rownames(limma_wholeProt_input_raw) <- proteome_Th17vTreg_base$prot.names
limma_wholeProt_input_raw <- limma_wholeProt_input_raw[proteome_Th17vTreg_base$max_pep >= max_pep_thresh_prot,]

experimental_design_prot <- data.frame(row.names=experiments_wholeProt,"Type"=str_extract(experiments_wholeProt, "WT|Bio"))
experimental_design_prot["condition"] <- factor(celltypes_prot, levels = c("Th17", "Treg"))
experimental_design_prot <- experimental_design_prot %>% rownames_to_column("label")
experimental_design_prot["technical_replicate"] <- str_extract(experiments_wholeProt, "[1-3]$")
experimental_design_prot["Exp"] <- str_extract(experiments_wholeProt, "Exp..[0-9][0-9]")
experimental_design_prot <- experimental_design_prot %>% mutate(replicate=paste(Exp, technical_replicate, sep="_"))
# summarized experiment object
bool_col <- which(colnames(proteome_Th17vTreg_base) %in% experiments_wholeProt)
# unique_prot <- make_unique(proteome_Th17vTreg_base, "prot.names", "Protein.Group", delim = ".")
proteome_Th17vTreg_base$name <- proteome_Th17vTreg_base$prot.names
proteome_Th17vTreg_base$ID <- proteome_Th17vTreg_base$Protein.Group
dep_input_raw <- make_se(proteome_Th17vTreg_base[proteome_Th17vTreg_base$max_pep >= max_pep_thresh,], bool_col, experimental_design_prot)

# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(dep_input_raw) %>%
  group_by(name, condition, Exp) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(dep_input_raw) %in% proteins_MNAR

# Perform a mixed imputation
dep_input_raw_imp <- impute(
  dep_input_raw, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "min") # imputation function for MNAR


limma_wholeProt_input <- assay(dep_input_raw_imp)
colnames(limma_wholeProt_input) <- colnames(limma_wholeProt_input_raw)
rownames(limma_wholeProt_input) <- rownames(limma_wholeProt_input_raw)


# design for limma
Cell <- factor(anno_wholeProt$Cell)
design_prot <- model.matrix(~0+Cell, data = anno_wholeProt, level=Th17)
contr.matrix <- makeContrasts(
  Th17vsTreg = CellTh17 - CellTreg,
  levels = design_prot)
# contr.matrix
# fit model to data
fit <- lmFit(limma_wholeProt_input, design_prot) #, robust=T, trend=T)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = T, robust = T)

# # qc fit
# plotMD(limma_wholeProt_input_raw)
# abline(0,0)
# plotMD(limma_wholeProt_input)
# abline(0,0)
# plotMD(fit)
# abline(0,0)
# plotMDS(limma_wholeProt_input)
# abline(0,0)
# voom(limma_wholeProt_input, design = design_prot, plot = T)
# abline(0,0)
# # check ebayes line
# qqt(fit$t,df=fit$df.prior+fit$df.residual,pch=16,cex=0.2)
# abline(0,1)


# volcano plot
results_prot <- as.data.frame(topTable(fit, number = dim(limma_wholeProt_input_raw)[1]))
results_prot_comp <- as.data.frame(topTable(fit, number = dim(limma_wholeProt_input_raw)[1]))
# custom colors for Treg and Th17 in volcano plot
keyvals <- ifelse(
  results_prot_comp$logFC < -logFC_prot & results_prot_comp$adj.P.Val < p_val_thresh, '#3E8EB9', # treg
  ifelse(results_prot_comp$logFC > logFC_prot & results_prot_comp$adj.P.Val < p_val_thresh, '#9C0E0F', # th17
         'darkgrey'))
keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == '#3E8EB9'] <- 'Treg_enriched'
names(keyvals)[keyvals == '#9C0E0F'] <- 'Th17_enriched'
names(keyvals)[keyvals == 'darkgrey'] <- 'ns'

EnhancedVolcano(results_prot_comp,
                lab = rownames(results_prot_comp),
                x = "logFC", 
                y = "adj.P.Val",
                selectLab = rownames(results_prot_comp)[which(names(keyvals) %in% c('Treg_enriched', 'Th17_enriched'))],
                colCustom = keyvals,
                pCutoff = p_val_thresh,
                FCcutoff = logFC_prot,
                pointSize = 1.5,
                labSize = 4.0,
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                title = "Proteome limma results",
                subtitle = "imputed data (min + knn)") 

ggsave(file = file.path(outdir, '/volcanoplot_wholeprot.pdf'), width=2000, height=2000, units = "px", dpi=300)

# helper functions for counting nas in dfs
count_na <- function(x) sum(is.na(x) | is.infinite(x))    # helper function
count_not_na <- function(x) sum(!is.na(x))    # helper function

th17_runfilter <- c(limma_wholeProt_input_raw[anno_wholeProt$Cell == "Th17"] %>% apply(1, count_na)) > max_NAs
treg_runfilter <- c(limma_wholeProt_input_raw[anno_wholeProt$Cell == "Treg"] %>% apply(1, count_na)) > max_NAs

results_prot <- results_prot[rownames(limma_wholeProt_input_raw[!(th17_runfilter & treg_runfilter),]), ] # remove Proteins, that do not have at least one condition with more than 8 non-NA samples 
results_prot <- results_prot[results_prot$adj.P.Val < p_val_thresh,] # p_val_thresh? more stringent like in original ttest
results_prot <- results_prot[abs(results_prot$logFC) >= logFC_prot,] # pulldown higher than 0.5! mind. 1? 
write.xlsx(results_prot, file.path(outdir, "/results_limma_proteome.xlsx"), rowNames=T) # prot names are missing :S



################################################################################
################################################################################
#                                Interactome                                   #      
################################################################################
################################################################################
# missing data mixed imp for Pulldown limma results with DEP imputation.
# Interactome is filtered for:
# maxpep >= 3 
# runfilter: 
# BIO_runfilter_treg <- c(limma_Treg_input_raw[anno_Treg$Type == "Bio"] %>% apply(1, count_na)) > 1
# results_paired_treg <- results_paired_treg[rownames(limma_Treg_input_raw[!(BIO_runfilter_treg),]), ] # remove Proteins, that do not have at least one condition with more than 8 non-NA samples 
# this means that we generate a boolean vector for the celltype of interest and search only in the Bio cells:
# holding a TRUE evertime a Protein has more than 1 NA in either of the groups.
# The vector is negated resulting in a dataset
# with only Proteins holdiung more than 8 non-NA values in the raw intensities of Bio.
# adj.P <= 0.05
# logFC > 1

folder.iact <- file.path("newest_data", "time_series", "Interaktom")
file.iact <- "2022-144_Pulldowns_Idefix_MaxLFQ_Intensities_AG.xlsx"

pulldown_combi <- read.xlsx(file.path(folder.iact, file.iact), sheet = 4, colNames = F)
names(pulldown_combi) <- pulldown_combi[2,]
names(pulldown_combi) <- make.names(names(pulldown_combi), unique = T)

pulldown_combi <- pulldown_combi %>% filter(Protein.Group != "Protein.Group")
pulldown_combi$prot.names <-  str_remove_all(pulldown_combi$prot.names, "_MOUSE")
# pulldown_combi <- pulldown_combi %>% dplyr::select(-sparklines) %>% filter(`SIGNIFICANT_BH.+.Log2.Filter` != "common")
pulldown_combi_base <- read.xlsx(file.path(folder.iact, file.iact), sheet = 2)
pulldown_combi_base <- pulldown_combi_base %>% dplyr::select(-X1)
pulldown_combi_base$prot.names <-  str_remove_all(pulldown_combi_base$prot.names, "_MOUSE")
# pulldown_combi_base[is.na(pulldown_combi_base)] <- 0
experiments_combi <- colnames(pulldown_combi)
experiments_combi <- experiments_combi[grepl("Exp", experiments_combi)]
anno_combi <- data.frame(row.names=experiments_combi,"Type"=str_extract(experiments_combi, "WT|Bio"), "Exp"=str_extract(experiments_combi, "Exp..[0-9]*"))
celltypes_prot <- str_extract(experiments_combi, "Th17|Treg")
anno_combi["Cell"] <- factor(celltypes_prot, levels = c("Th17", "Treg"))
anno_combi_row <- data.frame(row.names =pulldown_combi$prot.names,"max_pep"=unlist(lapply(pulldown_combi$max_pep, as.numeric)))
# remove Exp..59 from the analysis
anno_combi <- anno_combi %>% filter(Exp != "Exp..59")

experiments_combi <- rownames(anno_combi)


experiments_Th17 <- rownames(anno_combi %>% filter(Cell == "Th17"))
experiments_Treg <- rownames(anno_combi %>% filter(Cell == "Treg"))
pulldown_Th17 <- pulldown_combi %>% dplyr::select(prot.names, Protein.Group, all_of(experiments_Th17), `SIGNIFICANT_BH...Log2.Filter.1`)
pulldown_Th17 <- pulldown_Th17 %>% filter(`SIGNIFICANT_BH...Log2.Filter.1` != "ns")
pulldown_Treg <- pulldown_combi %>% dplyr::select(prot.names, Protein.Group, all_of(experiments_Treg), `SIGNIFICANT_BH...Log2.Filter`)
pulldown_Treg <- pulldown_Treg %>% filter(`SIGNIFICANT_BH...Log2.Filter` != "ns")

pulldown_Th17_base <- pulldown_combi_base %>% dplyr::select(prot.names, Protein.Group, max_pep, all_of(experiments_Th17))
pulldown_Treg_base <- pulldown_combi_base %>% dplyr::select(prot.names, Protein.Group, max_pep, all_of(experiments_Treg))

anno_Treg <- anno_combi %>% filter(Cell == "Treg")
anno_Th17 <- anno_combi %>% filter(Cell == "Th17")


################################################################################
# limmma pulldown Treg
Type <- factor(anno_Treg$Type)
Exp <- factor(anno_Treg$Exp)
design_Treg <- model.matrix(~0+Type+Exp, data = anno_Treg, level=Bio)
contr.matrix <- makeContrasts(
  BiovsWT = TypeBio - TypeWT,
  levels = design_Treg)
contr.matrix

# make input_raw
limma_Treg_input_raw <- log2(pulldown_Treg_base[experiments_Treg])
rownames(limma_Treg_input_raw) <- pulldown_Treg_base$prot.names
limma_Treg_input_raw <- limma_Treg_input_raw[pulldown_Treg_base$max_pep >= max_pep_thresh,]

# the block below is neccesary to fit the data into DEP
experimental_design_Treg <- data.frame(row.names=experiments_Treg,"condition"=str_extract(experiments_Treg, "WT|Bio"))
experimental_design_Treg <- experimental_design_Treg %>% rownames_to_column("label")
experimental_design_Treg["technical_replicate"] <- str_extract(experiments_Treg, "[1-3]$")
experimental_design_Treg["Exp"] <- str_extract(experiments_Treg, "Exp..[0-9][0-9]")
experimental_design_Treg <- experimental_design_Treg %>% mutate(replicate=paste(Exp, technical_replicate, sep="_"))
# summarized experiment object
bool_col <- which(colnames(pulldown_Treg_base) %in% experiments_Treg)
unique_prot <- make_unique(pulldown_Treg_base, "prot.names", "Protein.Group", delim = ";")
dep_input_raw <- make_se(unique_prot[unique_prot$max_pep >= max_pep_thresh,], bool_col, experimental_design_Treg)

# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(dep_input_raw) %>%
  group_by(name, condition, Exp) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(dep_input_raw) %in% proteins_MNAR

# Perform a mixed imputation
dep_input_raw_imp <- impute(
  dep_input_raw, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "min") # imputation function for MNAR

limma_treg_input <- assay(dep_input_raw_imp)
colnames(limma_treg_input) <- colnames(limma_Treg_input_raw)
rownames(limma_treg_input) <- rownames(limma_Treg_input_raw)

# input_raw matrix is done =>
################################################################################
# fit the model
fit <- lmFit(limma_treg_input, design_Treg )#, robust=T, trend=T)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = T, robust = T)
# plotSA(fit)
# qc fit
# plotMD(limma_treg_input)
# abline(0,0)
# plotMDS(limma_treg_input)
# abline(0,0)
# voom(limma_treg_input, design = design_Treg, plot = T)
# abline(0,0)

results_paired_treg <- as.data.frame(topTable(fit, number = dim(limma_Treg_input_raw)[1]))
results_paired_treg_comp <- as.data.frame(topTable(fit, number = dim(limma_Treg_input_raw)[1]))
# custom colors for Treg and Th17 in volcano plot
keyvals <- ifelse(
  results_paired_treg_comp$logFC > logFC_iact & results_paired_treg_comp$adj.P.Val < p_val_thresh , '#3E8EB9', # treg
  'darkgrey')
keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == '#3E8EB9'] <- 'Treg_interactor'
names(keyvals)[keyvals == 'darkgrey'] <- 'ns'

EnhancedVolcano(results_paired_treg_comp,
                lab = rownames(results_paired_treg_comp),
                x = "logFC", 
                y = "adj.P.Val",
                selectLab = rownames(results_paired_treg_comp)[which(names(keyvals) %in% c('Treg_interactor'))],
                colCustom = keyvals,
                # outcomment line below for complete volcano
                xlim = c(0, max(results_paired_treg_comp[["logFC"]], na.rm = TRUE) +1.5),
                pCutoff = p_val_thresh,
                FCcutoff = logFC_iact,
                pointSize = 1.5,
                labSize = 4.0,
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                title = "Treg limma results",
                subtitle = "imputed data (min + knn)") 
ggsave(file = file.path(outdir,'/paired_volcanoplot_treg.pdf'), width=2000, height=2000, units = "px", dpi=300)


BIO_runfilter_treg <- c(limma_Treg_input_raw[anno_Treg$Type == "Bio"] %>% apply(1, count_na)) > max_NAs
# WT_runfilter <- c(limma_Treg_input_raw[anno_Treg$Type == "WT"] %>% apply(1, count_na)) > 1

results_paired_treg <- results_paired_treg[rownames(limma_Treg_input_raw[!(BIO_runfilter_treg),]), ] # remove Proteins, that do not have at least one condition with more than 8 non-NA samples 
results_paired_treg <- results_paired_treg[results_paired_treg$adj.P.Val < p_val_thresh,] # p_val_thresh? more stringent like in original ttest
results_paired_treg <- results_paired_treg[results_paired_treg$logFC >= logFC_iact,] # pulldown higher than 0.5! mind. 1? 
write.xlsx(results_paired_treg, file.path(outdir,"/paired_results_limma_iact_Treg.xlsx"), rowNames=T) # prot names are missing :S


################################################################################
##### pulldown Th17

# built design matrix
Type <- factor(anno_Th17$Type) # BIO vs WT
Exp <- factor(anno_Th17$Exp) # the paired part of the analysis
design_th17 <- model.matrix(~0+Type+Exp, data = anno_Th17, level=Bio)
contr.matrix <- makeContrasts(
  BiovsWT = TypeBio - TypeWT,
  levels = design_th17)
contr.matrix

# make input_raw
limma_th17_input_raw <- log2(pulldown_Th17_base[experiments_Th17])
rownames(limma_th17_input_raw) <- pulldown_Th17_base$prot.names
limma_th17_input_raw <- limma_th17_input_raw[pulldown_Th17_base$max_pep >= max_pep_thresh,]

# the block below is neccesary to fit the data into DEP
experimental_design_th17 <- data.frame(row.names=experiments_Th17,"condition"=str_extract(experiments_Th17, "WT|Bio"))
experimental_design_th17 <- experimental_design_th17 %>% rownames_to_column("label")
experimental_design_th17["technical_replicate"] <- str_extract(experiments_Th17, "[1-3]$")
experimental_design_th17["Exp"] <- str_extract(experiments_Th17, "Exp..[0-9][0-9]")
experimental_design_th17 <- experimental_design_th17 %>% mutate(replicate=paste(Exp, technical_replicate, sep="_"))
# summarized experiment object
bool_col <- which(colnames(pulldown_Th17_base) %in% experiments_Th17)
unique_prot <- make_unique(pulldown_Th17_base, "prot.names", "Protein.Group", delim = ";")
dep_input_raw <- make_se(unique_prot[unique_prot$max_pep >= max_pep_thresh,], bool_col, experimental_design_th17)

# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(dep_input_raw) %>%
  group_by(name, condition, Exp) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()
# Get a logical vector
MNAR <- names(dep_input_raw) %in% proteins_MNAR

# Perform a mixed imputation
dep_input_raw_imp <- impute(
  dep_input_raw, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "min") # imputation function for MNAR

limma_th17_input <- assay(dep_input_raw_imp)
colnames(limma_th17_input) <- colnames(limma_th17_input_raw)
rownames(limma_th17_input) <- rownames(limma_th17_input_raw)
# input_raw matrix is done =>
################################################################################
# fit the model
fit <- lmFit(limma_th17_input, design_th17 )#, robust=T, trend=T)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = T, robust = T)
# plotSA(fit)
# extract results
results_paired_th17 <- as.data.frame(topTable(fit, number = dim(limma_th17_input_raw)[1]))
results_paired_th17_comp <- as.data.frame(topTable(fit, number = dim(limma_th17_input_raw)[1]))
# custom colors for Treg and Th17 in volcano plot
keyvals <- ifelse(
  results_paired_th17_comp$logFC > logFC_iact & results_paired_th17_comp$adj.P.Val < p_val_thresh , '#9C0E0F', # treg
  'darkgrey')
keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == '#9C0E0F'] <- 'Th17_interactor'
names(keyvals)[keyvals == 'darkgrey'] <- 'ns'
# volcano plot
EnhancedVolcano(results_paired_th17_comp,
                lab = rownames(results_paired_th17_comp),
                x = "logFC", 
                y = "adj.P.Val",
                selectLab = rownames(results_paired_th17_comp)[which(names(keyvals) %in% c('Th17_interactor'))],
                colCustom = keyvals,
                # outcomment line below for complete volcano
                xlim = c(0, max(results_paired_th17_comp[["logFC"]], na.rm = TRUE) +1.5),
                pCutoff = p_val_thresh,
                FCcutoff = logFC_iact,
                pointSize = 1.5,
                labSize = 4.0,
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                title = "Th17 limma results",
                subtitle = "imputed data (min + knn)") 
ggsave(file = file.path(outdir,'/paired_volcanoplot_th17.pdf'), width=2000, height=2000, units = "px", dpi=300)

BIO_runfilter_th17 <- c(limma_th17_input_raw[anno_Th17$Type == "Bio"] %>% apply(1, count_na)) > max_NAs
# WT_runfilter <- c(limma_th17_input_raw[anno_Th17$Type == "WT"] %>% apply(1, count_na)) > 1

results_paired_th17 <- results_paired_th17[rownames(limma_th17_input_raw[!(BIO_runfilter_th17),]), ] # remove Proteins, that do not have at least one condition with more than 7 nonNA samples 
results_paired_th17 <- results_paired_th17[results_paired_th17$adj.P.Val < p_val_thresh,] # p_val_thresh? more stringent like in original ttest
results_paired_th17 <- results_paired_th17[results_paired_th17$logFC >= logFC_iact,] # pulldown higher than 0.5! mind. 1? 
write.xlsx(results_paired_th17, file.path(outdir,"/paired_results_limma_iact_th17.xlsx"), row.names=T) # prot names are missing :S


################################################################################
# combi
limma_combi_input_raw <- log2(pulldown_combi_base[experiments_combi])
rownames(limma_combi_input_raw) <- pulldown_combi_base$prot.names
limma_combi_input_raw <- limma_combi_input_raw[pulldown_combi_base$max_pep >= max_pep_thresh,]

experimental_design_combi <- data.frame(row.names=experiments_combi,"condition"=str_extract(experiments_combi, "WT|Bio"))
experimental_design_combi <- experimental_design_combi %>% rownames_to_column("label")
experimental_design_combi["technical_replicate"] <- str_extract(experiments_combi, "[1-3]$")
experimental_design_combi["Exp"] <- str_extract(experiments_combi, "Exp..[0-9][0-9]")
experimental_design_combi["Cell"] <- celltypes_prot
experimental_design_combi <- experimental_design_combi %>% mutate(replicate=paste(Exp, technical_replicate, sep="_"))
# summarized experiment object
bool_col <- which(colnames(pulldown_combi_base) %in% experiments_combi)
unique_prot <- make_unique(pulldown_combi_base, "prot.names", "Protein.Group", delim = ";")
dep_input_raw <- make_se(unique_prot[unique_prot$max_pep >= max_pep_thresh,], bool_col, experimental_design_combi)

# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(dep_input_raw) %>%
  group_by(name, condition, Exp) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(dep_input_raw) %in% proteins_MNAR

# Perform a mixed imputation
dep_input_raw_imp <- impute(
  dep_input_raw, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "min") # imputation function for MNAR

limma_combi_input <- assay(dep_input_raw_imp)
colnames(limma_combi_input) <- colnames(limma_combi_input_raw)
rownames(limma_combi_input) <- rownames(limma_combi_input_raw)

Type <- factor(anno_combi$Type)
Cell <- factor(anno_combi$Cell)
Exp <- factor(anno_combi$Exp)
design_combi <- model.matrix(~0+Type+Cell+Exp, data = anno_combi, level=Bio)
contr.matrix <- makeContrasts(
  BiovsWT = TypeBio - TypeWT,
  levels = design_combi)
contr.matrix

# fit the model
fit <- lmFit(limma_combi_input, design_combi )#, robust=T, trend=T)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = T, robust = T)
# plotSA(fit)

results_paired_combi <- as.data.frame(topTable(fit, number = dim(limma_combi_input)[1]))
results_paired_combi_comp <- as.data.frame(topTable(fit, number = dim(limma_combi_input)[1]))
# volcano plot
EnhancedVolcano(results_paired_combi_comp,
                lab = rownames(results_paired_combi_comp),
                x = "logFC", 
                y = "adj.P.Val",
                # outcomment line below for complete volcano
                xlim = c(0, max(results_paired_combi_comp[["logFC"]], na.rm = TRUE) +1.5),
                pCutoff = p_val_thresh,
                FCcutoff = logFC_iact,
                pointSize = 1.5,
                labSize = 4.0,
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                title = "Combined limma results",
                subtitle = "imputed data (min + knn), Bio vs WT,\nwith Th17 and Treg as confounding factor") + 
  scale_colour_manual(values=cbPalette)
ggsave(file = file.path(outdir,'/paired_volcanoplot_combi.pdf'), width=2000, height=2000, units = "px", dpi=300)

# results_paired_combi <- results_paired_combi[rownames(limma_combi_input_raw[!(BIO_runfilter & WT_runfilter),]), ] # remove Proteins, that do not have at least one condition with more than 8 non-NA samples 
results_paired_combi <- results_paired_combi[results_paired_combi$adj.P.Val < p_val_thresh,] # p_val_thresh? more stringent like in original ttest
results_paired_combi <- results_paired_combi[results_paired_combi$logFC >= logFC_iact,] # pulldown higher than 0.5! mind. 1? 
write.xlsx(results_paired_combi, file.path(outdir,"/paired_results_limma_iact_combi.xlsx"), rowNames=T) # prot names are missing :S


################################################################################
#### output
data.folder <- "newest_data"

out_excel <- read.xlsx(file.path(data.folder, "2021-143_MaxLFQ_Intensities_Th17_TReg_Pulldown_together_MBR_nocrossnorm_Core4_UD.xlsx"), sheet = 3, colNames = F)
names(out_excel) <- out_excel[2,]
names(out_excel) <- make.names(names(out_excel), unique = T)
out_excel$prot.names <- str_remove_all(out_excel$prot.names, "_MOUSE")

th17_out <- as_tibble(results_paired_th17_comp, rownames = "prot.names") %>% 
  dplyr::select("prot.names", "logFC", "adj.P.Val") %>% 
  mutate(interactor = ifelse(prot.names %in% rownames(results_paired_th17), "th17_interactor", "ns"))
treg_out <- as_tibble(results_paired_treg_comp, rownames = "prot.names") %>% 
  dplyr::select("prot.names", "logFC", "adj.P.Val") %>% 
  mutate(interactor = ifelse(prot.names %in% rownames(results_paired_treg), "treg_interactor", "ns"))

out_excel <- out_excel %>% left_join(th17_out, by="prot.names") %>% 
  left_join(treg_out, by="prot.names", suffix=c("_Th17", "_Treg")) %>%
  left_join(as_tibble(BIO_runfilter_th17, rownames="prot.names") %>% dplyr::rename(th17_Bio_runfilter_NA=value)) %>%
  left_join(as_tibble(BIO_runfilter_treg, rownames="prot.names") %>% dplyr::rename(treg_Bio_runfilter_NA=value))

out_excel[2, c("logFC_Th17", "adj.P.Val_Th17", "interactor_Th17", "logFC_Treg", "adj.P.Val_Treg", "interactor_Treg", "th17_Bio_more_than_1_NA", "treg_Bio_more_than_1_NA")] <- 
  c("logFC_Th17", "adj.P.Val_Th17", "interactor_Th17", "logFC_Treg", "adj.P.Val_Treg", "interactor_Treg", "th17_Bio_more_than_1_NA", "treg_Bio_more_than_1_NA")
write.xlsx(out_excel, file.path(outdir,"/paired_mega_table_iact.xlsx"), rowNames=F, colNames=F)

# combi

out_excel <- read.xlsx(file.path(data.folder, "2021-143_MaxLFQ_Intensities_Th17_TReg_Pulldown_together_MBR_nocrossnorm_Core4_UD.xlsx"), sheet = 3, colNames = F)
names(out_excel) <- out_excel[2,]
names(out_excel) <- make.names(names(out_excel), unique = T)
out_excel$prot.names <- str_remove_all(out_excel$prot.names, "_MOUSE")

combi_out <- as_tibble(results_paired_combi_comp, rownames = "prot.names") %>% 
  dplyr::select("prot.names", "logFC", "adj.P.Val") %>% 
  mutate(interactor = ifelse(prot.names %in% rownames(results_paired_combi), "combi_interactor", "ns"))

colnames(combi_out) <-  c("prot.names", "logFC_combi", "adj.P.Val_combi", "interactor_combi")

out_excel <- out_excel %>% left_join(th17_out %>% 
                                       left_join(treg_out, by="prot.names", suffix=c("_Th17", "_Treg")) %>%
                                       left_join(combi_out, by="prot.names"),
                                     by="prot.names") %>%
  left_join(as_tibble(BIO_runfilter_th17, rownames="prot.names") %>% dplyr::rename(th17_Bio_more_than_1_NA=value)) %>%
  left_join(as_tibble(BIO_runfilter_treg, rownames="prot.names") %>% dplyr::rename(treg_Bio_more_than_1_NA=value))

out_excel[2, c("logFC_Th17", "adj.P.Val_Th17", "interactor_Th17", "logFC_Treg", "adj.P.Val_Treg", "interactor_Treg", "logFC_combi", "adj.P.Val_combi", "interactor_combi", "th17_Bio_more_than_1_NA", "treg_Bio_more_than_1_NA")] <- 
  c("logFC_Th17", "adj.P.Val_Th17", "interactor_Th17", "logFC_Treg", "adj.P.Val_Treg", "interactor_Treg", "logFC_combi", "adj.P.Val_combi", "interactor_combi", "th17_Bio_more_than_1_NA", "treg_Bio_more_than_1_NA")

write.xlsx(out_excel, file.path(outdir,"/paired_mega_table_iact_with_combi.xlsx"), rowNames=F, colNames=F)


# write imputed data to file
write.xlsx(as.data.frame(limma_th17_input), file.path(outdir,"/th17_data_imputed.xlsx"), rowNames=T)
write.xlsx(as.data.frame(limma_treg_input), file.path(outdir,"/treg_data_imputed.xlsx"), rowNames=T)


################################################################################
################################################################################
#                                Interactome:                                  #      
#                           th17 vs treg bio only:                             #
################################################################################
################################################################################

experiments_combi_bio <- experiments_combi[grepl("Bio", experiments_combi)]
subset_bio <- as.data.frame(limma_combi_input)[experiments_combi_bio]
# get only prots that were significant in either th17 or Treg between WT and Bio
subset_bio <- subset_bio %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% rownames(results_paired_th17) | prot.names %in% rownames(results_paired_treg)) %>%
  column_to_rownames("prot.names")

anno_subset <- anno_combi[experiments_combi_bio,]

Cell <- factor(anno_subset$Cell)
Exp <- factor(anno_subset$Exp)
#design_subset <- model.matrix(~0+Cell+Exp, data = anno_subset, level=Bio)
design_subset <- model.matrix(~0+Cell+Exp, data = anno_subset, level=Bio)
contr.matrix <- makeContrasts(
  Th17vsTreg = CellTh17 - CellTreg,
  levels = design_subset)
contr.matrix

# fit the model
fit <- lmFit(subset_bio, design_subset )#, robust=T, trend=T)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = T, robust = T)
# plotSA(fit)

# qc fit
# plotMD(subset_bio)
# abline(0,0)
# plotMDS(subset_bio)
# abline(0,0)
# voom(subset_bio, design = design_subset, plot = T)
# abline(0,0)


results_paired_subset <- as.data.frame(topTable(fit, number = dim(subset_bio)[1]))
results_paired_subset_comp <- as.data.frame(topTable(fit, number = dim(subset_bio)[1]))

# custom colors for Treg and Th17 in volcano plot
keyvals <- ifelse(
  results_paired_subset_comp$logFC < -logFC_prot & results_paired_subset_comp$adj.P.Val < p_val_thresh, '#3E8EB9', # treg
  ifelse(results_paired_subset_comp$logFC > logFC_prot & results_paired_subset_comp$adj.P.Val < p_val_thresh, '#9C0E0F', # th17
         'darkgrey'))
keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == '#3E8EB9'] <- 'Treg_enriched'
names(keyvals)[keyvals == '#9C0E0F'] <- 'Th17_enriched'
names(keyvals)[keyvals == 'darkgrey'] <- 'ns'
# volcano plot
EnhancedVolcano(results_paired_subset_comp,
                lab = rownames(results_paired_subset_comp),
                x = "logFC", 
                y = "adj.P.Val",
                selectLab = rownames(results_paired_subset_comp)[which(names(keyvals) %in% c('Treg_enriched', 'Th17_enriched'))],
                colCustom = keyvals,
                # outcomment line below for complete volcano
                # xlim = c(0, max(results_paired_subset_comp[["logFC"]], na.rm = TRUE) +1.5),
                ylim = c(0, max(-log10(results_paired_subset_comp[["adj.P.Val"]]), na.rm = TRUE) + 0.5),
                pCutoff = p_val_thresh,
                FCcutoff = logFC_prot,
                pointSize = 1.5,
                labSize = 4.0,
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                title = "Results Pulldown Bio Th17 vs Treg",
                subtitle = "subset of the Proteins used that were\nsignificant in either Treg or Th17.",
                # drawConnectors = TRUE,
                # widthConnectors = 0.75
)
ggsave(file = file.path(outdir,'/paired_volcanoplot_bioonly.pdf'), width=3000, height=3000, units = "px", dpi=300)

# results_paired_subset <- results_paired_subset[results_paired_subset$adj.P.Val < p_val_thresh,] # p_val_thresh? more stringent like in original ttest
# results_paired_subset <- results_paired_subset[results_paired_subset$logFC >= 1,] # pulldown higher than 0.5! mind. 1? 
write.xlsx(results_paired_subset %>%
             rownames_to_column("prot.names") %>%
             left_join(pulldown_combi_base %>% dplyr::select("prot.names", "Protein.Group")),
           file.path(outdir,"/paired_results_limma_iact_bioonly.xlsx"))


################################################################################
################################################################################
#                               Visualization                                  #      
################################################################################
################################################################################

# This block produces heatmaps, PCAs and the logRatio plot as well as a Venn Diagram. 
# There is already a plethora of plots arouynd here, however, they can easily be
# adapted or other made similarly by simple changeing/adding the parameters/ filters
# or pluffing in new data.
# see:
# ggplot2: https://ggplot2.tidyverse.org/reference/ggplot.html
# pheatmap: https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf , https://slowkow.com/notes/pheatmap-tutorial/
# EnhancedVolcano: https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html


# set paleddo
cbPalette <- c( "#E69F00", "#56B4E9","#999999", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#D4AC0D", "#7D3C98")
# hm palette
mat_colors <- inferno(10)

# column color annotations!
colors_matrix <- list(
  Cell = c("#9C0E0F", "#3E8EB9"),
  Type = c("#56B4E9", "#E69F00")
)
names(colors_matrix$Cell) <- c("Th17", "Treg")
names(colors_matrix$Type) <- c("WT", "Bio")

anno_Th17_hm <-  subset(anno_Th17, select = -Exp)
anno_Treg_hm <-  subset(anno_Treg, select = -Exp)
anno_subset_hm <-  subset(anno_subset, select = -Exp)
anno_combi_hm <- subset(anno_combi, select = -Exp)
# dev.off() # this needs to stay here. I suspect that ggsave and dev.off have some problems with each other, as the volcano plopts above are plotted correctly. 
# if removed heatmap_interactome_Th17.pdf is saved but with the logratio plot, for whatever reason.
# this is only apllicable for running the script via the console.
# if you run it interactively, feel free to drop it.


hm_inp <- as.matrix(limma_th17_input[rownames(results_paired_th17),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F, color=mat_colors,
               cluster_cols=T, annotation_col = anno_Th17_hm ,
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome_Th17.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome_Th17.csv"))


zscores <-  t(scale(t(limma_th17_input)))

hm_inp <- as.matrix(zscores[rownames(results_paired_th17),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F, color=mat_colors,
               cluster_cols=T, annotation_col = anno_Th17_hm ,
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome_Th17_zscore.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome_Th17_zscore.csv"))


hm_inp <- as.matrix(limma_treg_input[rownames(results_paired_treg),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T, annotation_col = anno_Treg_hm ,
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome_Treg.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome_Treg.csv"))

zscores <-  t(scale(t(limma_treg_input)))

hm_inp <- as.matrix(zscores[rownames(results_paired_treg),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T, annotation_col = anno_Treg_hm , 
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome_Treg_zscore.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome_Treg_zscore.csv"))

hm_inp <- as.matrix(subset_bio[rownames(results_paired_subset),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T, annotation_col = anno_subset_hm ,
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome_Th17vsTreg_subset.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome_Th17vsTreg_subset.csv"))

# ... as above with imputation
hm_inp <- as.matrix(limma_wholeProt_input[which(rownames(limma_wholeProt_input) %in% rownames(results_prot)),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T,
               annotation_col = anno_wholeProt,
               annotation_colors = colors_matrix,
               annotation_row = anno_wholeProt_row , annotation_legend = T)
ggsave(file = file.path(outdir,'/heatmap_Proteome.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)
# print(hm)
rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_Proteome.csv"))

# just the smallest of each
# plot gene tree
# pdf(file = file.path(outdir,'/heatmap_Proteome_tree.pdf'), width=1650, height=3000, res=300)
# plot(hm$tree_row)
# abline(h=40)
# dev.off()

# extract clusters
# clusts <- sort(cutree(hm$tree_row, h=40))
# #subset values
# hm_small <- hm_inp[c(names(clusts[which(clusts == 3)]),names( clusts[which(clusts == 4)]))  , ]
# hm <- pheatmap(hm_small, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
#                cluster_cols=T, annotation_col = anno_wholeProt,
#                annotation_colors = colors_matrix,
#                annotation_row = anno_wholeProt_row , annotation_legend = T)
# # print(hm)
# ggsave(file = file.path(outdir,'/heatmap_Proteome_sub.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)
# 
# rows_txt <- rownames(hm_small[hm$tree_row$order,])
# write.csv(rows_txt, file.path(outdir,"/heatmap_Proteome_sub.csv"))

## proteome zscore
# zscores <- limma_wholeProt_input %>% as.data.frame() %>% t() %>% mutate(across(everything(), ~ scale(.x))) %>% t()
zscores <-  t(scale(t(limma_wholeProt_input)))
hm_inp <- as.matrix(zscores[which(rownames(limma_wholeProt_input) %in% rownames(results_prot)),])
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T, annotation_col = anno_wholeProt,
               annotation_colors = colors_matrix,
               annotation_row = anno_wholeProt_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_Proteome_zscore.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_Proteome_zscore.csv"))

## HM of the complete interactome with imputed data
hm_inp <- as.matrix(limma_combi_input)
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T, annotation_col = anno_combi_hm, 
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome_complete.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome_complete.csv"))

## as above but filtered for significant genes
hm_inp <- as.matrix(limma_combi_input)[which(rownames(limma_combi_input) %in% c(rownames(results_paired_combi))) ,]
hm <- pheatmap(hm_inp, cluster_rows=T, show_rownames=F,show_colnames=F,annotation_names_row=F,color=mat_colors,
               cluster_cols=T, annotation_col = anno_combi_hm,
               annotation_colors = colors_matrix,
               annotation_row = anno_combi_row , annotation_legend = T)
# print(hm)
ggsave(file = file.path(outdir,'/heatmap_interactome.pdf'), plot=hm, width=1650, height=3000, units="px", dpi=300)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, file.path(outdir,"/heatmap_interactome.csv"))



################################################################################
# PCA all samples
# build combined data of all imputed data
prots_list <- unique(c(rownames(results_paired_th17),
                       rownames(results_paired_treg),
                       rownames(results_prot)))

all_omics_imp <- limma_wholeProt_input %>% as_tibble(rownames = "prot.names")
all_omics_imp <- all_omics_imp %>% 
  full_join(limma_th17_input %>% as_tibble(rownames = "prot.names"), by="prot.names",  suffix = c(".prot", ".iact")) %>%
  full_join(limma_treg_input %>% as_tibble(rownames = "prot.names"), by="prot.names",  suffix = c(".prot", ".iact")) %>%
  filter(prot.names %in% prots_list)
# annotation for all
experiments_all <- colnames(all_omics_imp)
experiments_all <- experiments_all[grepl("Exp", experiments_all)]
anno_omics_imp <- data.frame(row.names=experiments_all,"Type"=str_extract(experiments_all, "WT|Bio"))
anno_omics_imp["Experiment"] <- c(rep("Whole Proteome", times=length(experiments_wholeProt)) ,
                                  rep("Interactome", times=length(experiments_Th17)) ,
                                  rep("Interactome", times=length(experiments_Treg)) 
)
anno_omics_imp["Cell Type"] <- c(rep("Treg", times=length(experiments_wholeProt)/2) ,
                                 rep("Th17", times=length(experiments_wholeProt)/2) ,
                                 rep("Th17", times=length(experiments_Th17)) ,
                                 rep("Treg", times=length(experiments_Treg))
)

# get vector with info from which resp. experimenbt differential analysis stems
diff <- c()
iact_Th17 <- c()
iact_Treg <- c()
Th17_v_Treg <- c()

for (i in 1:length(all_omics_imp$prot.names)) {
  info <- ""
  if (all_omics_imp$prot.names[[i]] %in% rownames(results_prot)) {
    info <- paste0(info, " Th17 vs Treg")
    Th17_v_Treg <- c(Th17_v_Treg , "Th17 v Treg")
  } else {Th17_v_Treg <- c(Th17_v_Treg , "")}
  
  if (all_omics_imp$prot.names[[i]] %in% rownames(results_paired_th17)) {
    info <- paste0(info, " Th17 iact")
    iact_Th17 <- c(iact_Th17 , "Th17 iact")
  } else {iact_Th17 <- c(iact_Th17 , "")}
  if (all_omics_imp$prot.names[[i]] %in% rownames(results_paired_treg)) {
    info <- paste0(info, " Treg iact")
    iact_Treg <- c(iact_Treg , "Treg iact")
  } else {iact_Treg <- c(iact_Treg , "")}
  
  diff <- c(diff, info)
}
# names(diff) <- all_omics_imp$prot.names
anno_omics_imp_row <- data.frame(row.names =all_omics_imp$prot.names,
                                 # "max_pep"=all_omics_imp$max_pep, 
                                 # "max_PG.Q.val.prot"=all_omics_imp$`max_PG.Q.val.prot`,
                                 #"Differential"=diff
                                 "Proteome"=Th17_v_Treg,
                                 "Iact Th17"=iact_Th17,
                                 "iact Treg"=iact_Treg
)
anno_omics_imp_row_complete <- tibble("prot.names" =all_omics_imp$prot.names,
                                      #"max_pep"=all_omics_imp$max_pep, 
                                      # "max_PG.Q.val.prot"=all_omics_imp$`max_PG.Q.val.prot`,
                                      "Differential"=diff,
                                      "Proteome"=Th17_v_Treg,
                                      "Iact Th17"=iact_Th17,
                                      "iact Treg"=iact_Treg
)
# plot the actual pca
pca.table <- as.matrix(all_omics_imp %>% column_to_rownames("prot.names"))
res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)

pca_plot_df_all <-  as_tibble(res.pca$var$coord, rownames="id") %>% 
  left_join(anno_omics_imp %>% rownames_to_column("id"), by="id") 

pca.plot <- pca_plot_df_all %>% 
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell Type`, shape=`Type`)) +
  stat_ellipse(level = 0.99, aes(color = Experiment)) +
  scale_color_manual(values=c( "#E69F00", "#9C0E0F","#3E8EB9",  "#56B4E9"))
# save
ggsave(file.path(outdir,"/pca_all.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

# PCA iact imp
pca.table <- cbind(limma_th17_input, limma_treg_input) 
res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)

pca_plot_df_all <-  as_tibble(res.pca$var$coord, rownames="id") %>% 
  left_join(anno_combi %>% rownames_to_column("id"), by="id") 

pca.plot <- pca_plot_df_all %>% 
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
# save
ggsave(file.path(outdir,"/pca_interactome.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

pca.plot <- pca_plot_df_all %>%  
  dplyr::filter(Type == "Bio") %>%
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
# save
ggsave(file.path(outdir,"/pca_interactome_bio_only.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)
################################################################################
# pca Proteome
# high and low qualitles only
pca.table <- limma_wholeProt_input
# filter for variance
#  First find the desired quantile breaks for the entire matrix
qt <- quantile( pca.table , probs = c(0.01,0.99) , na.rm=T)
#  Next get a logical vector of the rows that have any values outside these breaks
rows <- apply( pca.table , 1 , function(x) any( x < qt[1] | x > qt[2] ) )
rows[is.na(rows)] <- FALSE # rmeove nas
#  Subset on this vector
pca.table <- pca.table[ rows , ]

res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)

pca_plot_df_proteome <-  as_tibble(res.pca$var$coord, rownames="id") %>% 
  left_join(anno_wholeProt %>% rownames_to_column("id"), by="id") 

pca.plot <- pca_plot_df_proteome %>% 
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
ggsave(file.path(outdir,"/pca_proteome_imp_only_high+low_quantile.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

# whole Prot full pca
pca.table <- limma_wholeProt_input
res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)

pca_plot_df_proteome <-  as_tibble(res.pca$var$coord, rownames="id") %>% 
  left_join(anno_wholeProt %>% rownames_to_column("id"), by="id") 

pca.plot <- pca_plot_df_proteome %>% 
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  # geom_text_repel(aes(label=id)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
ggsave(file.path(outdir,"/pca_proteome_imp.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

pca.plot <- pca_plot_df_proteome %>% 
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  geom_text_repel(aes(label=id), size=1.5, segment.size=0.1) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
ggsave(file.path(outdir,"/pca_proteome_imp_annotated.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

################################################################################
# pca Interactome
pca.table <- 2^cbind(limma_th17_input, limma_treg_input) 
#  First find the desired quantile breaks for the entire matrix
qt <- quantile( pca.table , probs = c(0.01,0.99) , na.rm=T)
#  Next get a logical vector of the rows that have any values outside these breaks
rows <- apply( pca.table , 1 , function(x) any( x < qt[1] | x > qt[2] ) )
rows[is.na(rows)] <- FALSE # rmeove nas
#  Subset on this vector
pca.table <- pca.table[ rows , ]

res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
pca_plot_df_interactome <-  as_tibble(res.pca$var$coord, rownames="id") %>% 
  left_join(anno_combi %>% rownames_to_column("id"), by="id")


pca.plot <- pca_plot_df_interactome %>% 
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9")) + 
  geom_text_repel(aes(label=id))
ggsave(file.path(outdir,"/pca_interactome_only_high+low_quantile.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

pca.plot <- pca_plot_df_interactome %>% 
  dplyr::filter(Exp != "Exp..26") %>%
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
ggsave(file.path(outdir,"/pca_interactome_only_high+low_quantile_noExp26.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)

pca.plot <- pca_plot_df_interactome %>% 
  dplyr::filter(Type == "Bio") %>%
  ggplot(aes(x= Dim.1, y= Dim.2)) +
  geom_point(aes(color = `Cell`, shape=`Type`)) +
  scale_color_manual(values=c("#9C0E0F","#3E8EB9"))
ggsave(file.path(outdir,"/pca_interactome_only_high+low_quantile_bio_only.pdf"), plot=pca.plot, width = 5, height = 5, dpi=300)


################################################################################
# venn diagram interactome vs interactome
wholeProt_prots <-  proteome_Th17vTreg_base %>%  
  dplyr::select("prot.names", all_of(experiments_wholeProt)) %>% 
  filter(prot.names %in%  rownames(results_prot)) %>% 
  pull("prot.names")

Th17_prots <-  pulldown_Th17_base %>%  
  dplyr::select("prot.names", all_of(experiments_Th17))  %>% 
  filter(prot.names %in%  rownames(results_paired_th17)) %>% 
  pull("prot.names") 

Treg_prots <-  pulldown_Treg_base %>%  
  dplyr::select("prot.names", all_of(experiments_Treg))  %>% 
  filter(prot.names %in%  rownames(results_paired_treg)) %>% 
  pull("prot.names") 

venn_input_raw <- list(
  Interactome_Th17=Th17_prots,
  Interactome_Treg=Treg_prots,
  Proteome_Th17vsTreg=wholeProt_prots
)

venno <- ggVennDiagram(venn_input_raw, edge_size=0) +
  scale_fill_gradient(low="steelblue",high = "darkorange") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(file.path(outdir,"/venn_interactome-vs-Proteome-ttest.pdf"), plot = venno, width = 7, height = 7, dpi = 300)

################################################################################
# log ratio plot
# proteome th17 vs treg against interactome th17 vs treg, bio only
# get log ratios => log2(mean(th17)/mean(treg))
# interactome
bio_th17 <- rownames(anno_Th17 %>% filter(Type == "Bio"))
bio_treg <- rownames(anno_Treg %>% filter(Type == "Bio"))
# no imp
# bioOnly_iact <- pulldown_Th17_base %>% dplyr::select("prot.names", all_of(bio_th17)) %>%
#   left_join(pulldown_Treg_base %>% dplyr::select("prot.names", all_of(bio_treg)), by="prot.names")
# imp
bioOnly_iact <- as_tibble(2^limma_th17_input, rownames="prot.names") %>% dplyr::select("prot.names", all_of(bio_th17)) %>%
  left_join(as_tibble(2^limma_treg_input , rownames="prot.names") %>% dplyr::select("prot.names", all_of(bio_treg)), by="prot.names")

bioOnly_iact_logratio <- log2(
  (rowMeans(bioOnly_iact %>% dplyr::select(all_of(bio_th17))) + 1) /
    (rowMeans(bioOnly_iact %>% dplyr::select(all_of(bio_treg))) + 1)
)
names(bioOnly_iact_logratio) <- bioOnly_iact$prot.names
# proteome
prot_th17 <- rownames( anno_omics_imp %>% filter(Experiment == "Whole Proteome" & `Cell Type`== "Th17"))
prot_th17 <- str_remove(prot_th17, ".prot")
prot_treg <- rownames( anno_omics_imp %>% filter(Experiment == "Whole Proteome" & `Cell Type`== "Treg"))
prot_treg <- str_remove(prot_treg, ".prot")
# no imp
# prot_logratio <- log2(
#   (rowMeans(proteome_Th17vTreg_base %>% dplyr::select(all_of(prot_th17))) + 1) /
#     (rowMeans(proteome_Th17vTreg_base %>% dplyr::select(all_of(prot_treg))) + 1)
# )
# names(prot_logratio) <- proteome_Th17vTreg_base$prot.names
# imp 
prot_logratio <- log2(
  (rowMeans(as_tibble(2^limma_wholeProt_input) %>% dplyr::select(all_of(prot_th17))) + 1) /
    (rowMeans(as_tibble(2^limma_wholeProt_input) %>% dplyr::select(all_of(prot_treg))) + 1)
)
names(prot_logratio) <- rownames(limma_wholeProt_input)

# filter only signififcant 
prots_list <- unique(c(rownames(results_paired_th17),
                       rownames(results_paired_treg),
                       rownames(results_prot)
)
)
bioOnly_iact_logratio <- bioOnly_iact_logratio[names(bioOnly_iact_logratio) %in% prots_list]
prot_logratio <- prot_logratio[names(prot_logratio) %in% prots_list]
### build df and plot
plot_df <- as.data.frame(bind_rows(Proteome=prot_logratio, Interactome=bioOnly_iact_logratio))
rownames(plot_df) <-c("log2ratio_Proteome", "log2ratio_Interactome")
plot_df <- as_tibble(t(plot_df), rownames="prot.names")

plot_df <- plot_df %>% left_join(anno_omics_imp_row_complete, by="prot.names")

plot_df <- as_tibble(plot_df %>% na.omit())

# plot_df <- plot_df %>% mutate("exceeds threshold"=log2ratio_Proteome >= 1 | 
#                                 log2ratio_Proteome<= -1 | 
#                                 log2ratio_Interactome >= 1 | 
#                                 log2ratio_Interactome <= -1)
plot_df <- plot_df %>% mutate("Interactome_only"= ifelse(grepl("vs", `Differential` ,), "No", "Yes"))

#### compute orthodgonal distances
# Compute the linear regression
fit <- lm(log2ratio_Interactome ~ log2ratio_Proteome, data=plot_df)

#finds endpoint for a perpendicular segment from the point (x0,y0) to the line
perp.segment.coord <- function(x0, y0, a = 0, b = 1) {
  x1 <- (x0 + b * y0 - a * b) / (1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
}

# Compute the endpoints of the orthogonal lines for each point
endpoints <- apply(cbind(plot_df$log2ratio_Proteome, plot_df$log2ratio_Interactome), 1, function(pt) {
  perp.segment.coord(pt[1], pt[2], fit$coefficients[1], fit$coefficients[2])
})
# Add the endpoints to the data frame
plot_df$xend <- sapply(endpoints, function(e) e$x1)
plot_df$yend <- sapply(endpoints, function(e) e$y1)
# compute euc dist
plot_df$distance <- sqrt((plot_df$log2ratio_Proteome - plot_df$xend)^2 + (plot_df$log2ratio_Interactome - plot_df$yend)^2)
# add threshold column
plot_df <- plot_df %>% mutate("exceeds threshold"=distance >= quantile(plot_df$distance, 0.9) ) # 0.9 quantile
#### add data column to plot_df that holds only prot names of proteins that exceed threshold.
plot_df <- plot_df %>%
  mutate(outside = ifelse(distance >= quantile(plot_df$distance, 0.9), prot.names, ""))

plot_df_quant <- plot_df[plot_df$distance > quantile(plot_df$distance, 0.9),]

plt <- plot_df %>% 
  ggplot(aes(x=log2ratio_Proteome, y=log2ratio_Interactome)) +
  geom_point( aes(color=`exceeds threshold`)) +
  geom_smooth(method = "lm", color="black") +
  geom_text_repel(aes(label=outside) )+ #, arrow = arrow(length = unit(0.02, "npc")) )+ #, hjust=1, vjust=1)+ 
  geom_segment(data = plot_df_quant, aes(x = log2ratio_Proteome, y = log2ratio_Interactome, xend = xend, yend = yend),
               linetype = "dashed", color = "blue") +
  # geom_hline(yintercept = 1) + 
  # geom_vline(xintercept = 1) +
  # geom_hline(yintercept = -1) + 
  # geom_vline(xintercept = -1) +
  scale_colour_manual(values=cbPalette) +
  guides(color=F)+
  theme_minimal() + coord_fixed()

ggsave(file.path(outdir,'/log2ratios-data-withProt.pdf'), plot = plt, width=12, height=10, dpi=300)
write.xlsx(plot_df, file.path(outdir,"/log2ratios-data-withProt.xlsx"))

# with limmas logFCs
bioOnly_iact_logratio <- results_paired_subset_comp$logFC
names(bioOnly_iact_logratio) <- rownames(results_paired_subset_comp)

prot_logratio <- results_prot_comp$logFC
names(prot_logratio) <- rownames(results_prot_comp)

# # filter only signififcant 
# prots_list <- unique(c(rownames(results_paired_th17),
#                        rownames(results_paired_treg),
#                        rownames(results_prot)
# )
# )
# bioOnly_iact_logratio <- bioOnly_iact_logratio[names(bioOnly_iact_logratio) %in% prots_list]
# prot_logratio <- prot_logratio[names(prot_logratio) %in% prots_list]
### build df and plot
plot_df <- as.data.frame(bind_rows(Proteome=prot_logratio, Interactome=bioOnly_iact_logratio))
rownames(plot_df) <-c("log2ratio_Proteome", "log2ratio_Interactome")
plot_df <- as_tibble(t(plot_df), rownames="prot.names")

plot_df <- plot_df %>% left_join(anno_omics_imp_row_complete, by="prot.names")

plot_df <- as_tibble(plot_df %>% na.omit())

# plot_df <- plot_df %>% mutate("exceeds threshold"=log2ratio_Proteome >= 1 | 
#                                 log2ratio_Proteome<= -1 | 
#                                 log2ratio_Interactome >= 1 | 
#                                 log2ratio_Interactome <= -1)
plot_df <- plot_df %>% mutate("Interactome_only"= ifelse(grepl("vs", `Differential` ,), "No", "Yes"))

#### compute orthodgonal distances
# Compute the linear regression
fit <- lm(log2ratio_Interactome ~ log2ratio_Proteome, data=plot_df)

#finds endpoint for a perpendicular segment from the point (x0,y0) to the line
perp.segment.coord <- function(x0, y0, a = 0, b = 1) {
  x1 <- (x0 + b * y0 - a * b) / (1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
}

# Compute the endpoints of the orthogonal lines for each point
endpoints <- apply(cbind(plot_df$log2ratio_Proteome, plot_df$log2ratio_Interactome), 1, function(pt) {
  perp.segment.coord(pt[1], pt[2], fit$coefficients[1], fit$coefficients[2])
})
# Add the endpoints to the data frame
plot_df$xend <- sapply(endpoints, function(e) e$x1)
plot_df$yend <- sapply(endpoints, function(e) e$y1)
# compute euc dist
plot_df$distance <- sqrt((plot_df$log2ratio_Proteome - plot_df$xend)^2 + (plot_df$log2ratio_Interactome - plot_df$yend)^2)
# add threshold column
plot_df <- plot_df %>% mutate("exceeds threshold"=distance >= quantile(plot_df$distance, 0.9) ) # 0.9 quantile
#### add data column to plot_df that holds only prot names of proteins that exceed threshold.
plot_df <- plot_df %>%
  mutate(outside = ifelse(distance >= quantile(plot_df$distance, 0.9), prot.names, ""))

plot_df_quant <- plot_df[plot_df$distance > quantile(plot_df$distance, 0.9),]

plt <- plot_df %>% 
  mutate(Interactome = ifelse(log2ratio_Interactome > 0.5, "Th17_enriched", ifelse(log2ratio_Interactome < -0.5, "Treg_enriched", "neither"))) %>% 
  ggplot(aes(x=log2ratio_Proteome, y=log2ratio_Interactome)) +
  # geom_point( aes(color=`exceeds threshold`)) +
  geom_point( aes(color=Interactome, shape=Proteome), size=3) +
  geom_smooth(method = "lm", color="black") +
  geom_text_repel(aes(label=outside))+ #, arrow = arrow(length = unit(0.02, "npc")) )+ #, hjust=1, vjust=1)+ 
  geom_segment(data = plot_df_quant, aes(x = log2ratio_Proteome, y = log2ratio_Interactome, xend = xend, yend = yend),
               linetype = "dashed", color = "black") +
  # geom_hline(yintercept = 1) + 
  # geom_vline(xintercept = 1) +
  # geom_hline(yintercept = -1) + 
  # geom_vline(xintercept = -1) +
  scale_colour_manual(values=list(neither="#999999", Th17_enriched="#9C0E0F", Treg_enriched="#2D789C")) +
  guides(color=F, shape=F)+
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  coord_fixed()


# 
ggsave(file.path(outdir,'/log2ratios-data-limma-iactanno.pdf'), plot = plt, width=12, height=10, dpi=300)
write.xlsx(plot_df, file.path(outdir,"/log2ratios-data-limma.xlsx"))