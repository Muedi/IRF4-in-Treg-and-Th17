######################################################################################################
# Script: KO-Data Analysis
# Author: Maximilian Sprang, Muedi
# Date: 20.09.2023
# Description: 
# This script uses Limma to analyse IRF4 KO Vs WT datasets
# We also compare day 0 vs day 3 to obtain information about Treg and Th17 differentiation.
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
library(sjmisc)
library(plotly)
library(pheatmap)
library(EnhancedVolcano)
library(stringr)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dendextend)


data.folder <- "newest_data"
# load meta_complete and dessolve conditions
meta_complete <- read_csv(file.path(data.folder, "2022-041_rawfile_mapping_Hauptmessung_IRF4WT_Ko.csv"))
meta_complete <- meta_complete %>% mutate(cell_diff = ifelse(str_detect(condition, "naive"), "naive", ifelse(str_detect(condition, "Th17"), "Th17", "Treg")))
meta_complete <- meta_complete %>% mutate(IRF4_stat = ifelse(str_detect(condition, "Ko"), "KO", "WT"))
meta_complete <- meta_complete %>% mutate(time_point =  ifelse(str_detect(condition, "D0"), "D0", "D3"))
meta_complete <- meta_complete %>% mutate(condi_real =  str_remove(condition, " n=[0-9]"))
meta_complete <- meta_complete %>% mutate(condi_real =  gsub(" ", "_", condi_real))
meta_complete <- as.data.frame(meta_complete)
rownames(meta_complete) <- meta_complete$identifier

condis <- unique(meta_complete$condi_real)

# load data
data_complete <- read_csv(file.path(data.folder, "2022-041_Hauptmessung_TTP_IRF4KovsWT_Intensities.csv"))
data_complete$prot.names <- str_remove(data_complete$prot.names, "_MOUSE")
# get experiments and non other columns
experiments_multi <- colnames(data_complete)
experiments_multi <- experiments_multi[grepl("IRF4", experiments_multi)]
non_exp_cols <- colnames(data_complete)
non_exp_cols <- non_exp_cols[!grepl("IRF4", non_exp_cols)]
# max pep filter
data_complete <- data_complete[data_complete$max_pep > 2,]
# remove rows with 2/3 of NA for this sample type
# data_complete = filter(data_complete, !apply( data_complete[experiments_multi], 1, function(x) sum(is.na(x)) >= dim(data_complete[samples])[2]/3*2 ) )

## Conditions
# "IRF4_Ko_1_naive_CD4+_D0" "IRF4_WT_1_naive_CD4+_D0"
# "IRF4_Ko_1_Th17_D3"       "IRF4_WT_1_Th17_D3"      
# "IRF4_Ko_1_TReg_D3"       "IRF4_WT_1_TReg_D3"   

### Comparisons:
naiveWT_vs_KO <- c("IRF4_WT_1_naive_CD4+_D0", "IRF4_Ko_1_naive_CD4+_D0")
Th17WT_vs_KO <- c("IRF4_WT_1_Th17_D3", "IRF4_Ko_1_Th17_D3")
TregWT_vs_KO <- c("IRF4_WT_1_TReg_D3", "IRF4_Ko_1_TReg_D3")
Th17KO_vs_TregWT <- c("IRF4_Ko_1_Th17_D3", "IRF4_WT_1_TReg_D3" )

# unpaired:
naiveWT_vs_Th17WT <- c("IRF4_WT_1_naive_CD4+_D0", "IRF4_WT_1_Th17_D3")
naiveWT_vs_Th17KO <- c("IRF4_WT_1_naive_CD4+_D0", "IRF4_Ko_1_Th17_D3")
naiveWT_vs_TregWT <- c("IRF4_WT_1_naive_CD4+_D0", "IRF4_WT_1_TReg_D3")
naiveWT_vs_TregKO <- c("IRF4_WT_1_naive_CD4+_D0", "IRF4_Ko_1_TReg_D3")


comparisons <- list(naiveWT_vs_KO,
                    Th17WT_vs_KO,
                    TregWT_vs_KO,
                    Th17KO_vs_TregWT,
                    naiveWT_vs_Th17WT,
                    naiveWT_vs_Th17KO,
                    naiveWT_vs_TregWT,
                    naiveWT_vs_TregKO)
names(comparisons) <- c("naiveWT_vs_KO",
                      "Th17WT_vs_KO",
                      "TregWT_vs_KO",
                      "Th17KO_vs_TregWT",
                      "naiveWT_vs_Th17WT",
                      "naiveWT_vs_Th17KO",
                      "naiveWT_vs_TregWT",
                      "naiveWT_vs_TregKO")

# outfiles
out <- data_complete


# palette fopr volcanos
cbPalette <- c( "#999999", "#56B4E9", "#009E73", "#E69F00")
# cbPalette <- inferno(5)

# save GO info
GO_res <- list()

### will be aloop in the future.
for (i in 1:length(comparisons)){
  samples <- meta_complete %>% filter((condi_real == comparisons[[i]][[1]]) |
                          (condi_real == comparisons[[i]][[2]])) %>%
                          pull(identifier)
  control <- str_split(names(comparisons[i]), "_vs_")[[1]][1]
  other <- str_split(names(comparisons[i]), "_vs_")[[1]][2]
  
  
  data <- data_complete %>% dplyr::select(all_of(non_exp_cols), all_of(samples)) 
  meta <- meta_complete %>% filter(identifier %in% samples)
  
  # filter rows, where in one Cell/genotype combimore than 3 NAs are present
  data_inp = data[samples] %>% as.data.frame()
  rownames(data_inp) <- data$prot.names
  data_inp <- data_inp[rownames(meta)]
  
  # run filter
  count_na <- function(x) sum(is.na(x))    # helper function
  count_not_na <- function(x) sum(!is.na(x))    # helper function

  condi1_more_then_9 <- c(data_inp[meta$condi_real == comparisons[[i]][[1]] ] %>% apply(1, count_na)) > 3
  condi2_more_then_9 <- c(data_inp[meta$condi_real == comparisons[[i]][[2]] ] %>% apply(1, count_na)) > 3
  
  combined_bool <- (condi1_more_then_9 & condi2_more_then_9)
  
  # filter 
  data <- data[!combined_bool, ]
  data_inp <- data_inp[!combined_bool, ]
  
  ##############################################################################
  ####################### mixed IMPUTATION #####################################
  # get summarized experiment
  experimental_design_mult <- data.frame(row.names=samples,"Type"=str_extract(samples, "WT|Ko"))
  experimental_design_mult["condition"] <- ifelse(grepl("D0", samples), "Naive",
                                                  ifelse(grepl("D3_TReg", samples), "Treg", "Th17"))
  experimental_design_mult <- experimental_design_mult %>% rownames_to_column("label")
  experimental_design_mult["Exp"] <- str_extract(samples, "[1-3]$")
  experimental_design_mult["technical_replicate"] <- str_sub(samples, 1, 9)
  experimental_design_mult <- experimental_design_mult %>% mutate(replicate=paste(Exp, technical_replicate, sep="_"))

  # summarized experiment object
  bool_col <- which(colnames(data) %in% samples)
  unique_mult <- make_unique(data, "prot.names", "Protein.Group", delim = ".")
  # the fauntion failed, but did not throw an error, so we will fill in name and id ourselves.
  unique_mult$name <- unique_mult$prot.names
  unique_mult$id <- unique_mult$protein.group
  dep_input <- make_se(unique_mult[unique_mult$max_pep >= 2,], bool_col, experimental_design_mult)

  # extract protein names with missing values
  # in all replicates of at least one condition
  proteins_MNAR <- get_df_long(dep_input) %>%
    group_by(name, condition, Exp) %>%
    summarize(nas = all(is.na(intensity))) %>%
    filter(nas) %>%
    pull(name) %>%
    unique()

  # get a logical vector
  MNAR <- names(dep_input) %in% proteins_MNAR

  # perform a mixed imputation
  dep_input_imp <- impute(
    dep_input,
    fun = "mixed",
    randna = !MNAR, # we have to define mar which is the opposite of MNAR
    mar = "knn", # imputation function for mar
    mnar = "min") # imputation function for MNAR


  limm_inp <- assay(dep_input_imp)
  colnames(limm_inp) <- samples

  # check distribution 
  # hist(limm_inp, breaks=100)
  # hist(as.matrix(log2(data_inp[samples])), breaks=100)
  
  limm_inp <- as.data.frame(limm_inp)
  
  # write out imputed data if wanted
  # write.xlsx(as.data.frame(limm_inp),
  #            paste0("output_new/mouse_knockout/imputed_data_",comparisons[[i]][[1]] ,"_vs_",comparisons[[i]][[2]] , ".xlsx"),
  #            rownames=t)
  
  
  
  ##############################################################################
  ############################### LIMMA 
  ##############################################################################
  
  # assign correct control
  cond <- relevel(cond, comparisons[[i]][[1]])
  design_base = model.matrix(~cond, data=limm_inp) # assign factors to the samples
  rownames(design_base) = colnames(limm_inp)
  colnames(design_base) = str_remove(colnames(design_base), "cond")
  colnames(design_base) = str_remove(colnames(design_base), "CD4[+]_")
  
  # cont = makeContrasts(
  #   naive_KO_WT=IRF4_Ko_1_naive_D0-IRF4_WT_1_naive_D0,
  # 
  #   levels=design_base)
  # 
  ## fit the model to the data:
  fit = lmFit(limm_inp, design_base)

  # contrast
  # fit = contrasts.fit(fit, cont)
  # # make statistics for the fit
  fit = eBayes(fit)
  # quick check the results for each contrast
  results <- decideTests(fit)
  summary(results)
  
  result_table = topTable(fit, adjust="BH", number=Inf)
  result_table <- result_table %>% rownames_to_column("prot.names")
  write.csv(result_table, paste0("output_new/mouse_knockout/results_", names(comparisons[i]), ".csv"),  row.names = T, col.names = T)
  
  # volcano
  plt <- EnhancedVolcano(result_table,
                  lab = result_table$prot.names,
                  x = "logFC", 
                  y = "adj.P.Val",
                  pCutoff = 0.01,
                  FCcutoff = 0.5,
                  pointSize = 1.5,
                  labSize = 4.0,
                  legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                  title = paste0(names(comparisons[i]), " results"),
                  subtitle = "imputed data (min + knn)") + 
    scale_colour_manual(values=cbPalette)
  
  ggsave(file = paste0('output_new/mouse_knockout/volcano_', names(comparisons[i]),'.pdf'),plot=plt, width=2000, height=2000, units = "px",  dpi=300)
  
  
  ##############################################################################
  ############################## GO enrichment 
  ##############################################################################

  protList <- result_table$t
  names(protList) <- data_complete[ match(result_table$prot.names, data_complete$prot.names),]$Protein.Group
  protList <- sort(protList, decreasing = T)
  
  ego <- gseGO(geneList      = protList,
                OrgDb        = org.Mm.eg.db,
                keyType      = "UNIPROT",
                ont          = "MF",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
  
  GO_res[names(comparisons[i])] = ego@result[1, "Description"]
  
  result_table["Protein.Group"] <- data_complete[ match(result_table$prot.names, data_complete$prot.names),]$Protein.Group
  result_table["GO_annos"] <- NA
  for (j in 1:3) {
      result_table <- result_table %>% 
      mutate("GO_annos" = ifelse(str_detect(ego@result[j,"core_enrichment"], Protein.Group), 
                                 ego@result[j,"Description"], 
                                 GO_annos))
  }
  
  # join to output
  result_table <- result_table %>% column_to_rownames("prot.names")
  result_table <- result_table %>% 
    dplyr::select("logFC", "adj.P.Val", "GO_annos") %>%
    mutate(signi = ifelse(adj.P.Val < 0.01 & logFC > 0.5, other, ifelse(
      adj.P.Val < 0.01 & logFC < -0.5, control, "ns"))
      )
  names(result_table) <- paste0(names(result_table), "_", names(comparisons[i]))
  result_table <- result_table %>% rownames_to_column("prot.names")
  out <- out %>% left_join(result_table, by="prot.names")
  
}

out <- out %>% mutate(Nas_per_prot = apply( out[experiments_multi], 1, function(x) sum(is.na(x)) )) 
out <- out %>% mutate(Nas_in_WT_naive = apply( out[rownames(meta[meta$condi_real == "IRF4_WT_1_naive_CD4+_D0",])], 1, function(x) sum(is.na(x)) )) 
out <- out %>% mutate(Nas_in_KO_naive = apply( out[rownames(meta[meta$condi_real == "IRF4_Ko_1_naive_CD4+_D0",])], 1, function(x) sum(is.na(x)) )) 
out <- out %>% mutate(Nas_in_WT_Th17 = apply( out[rownames(meta[meta$condi_real == "IRF4_WT_1_Th17_D3",])], 1, function(x) sum(is.na(x)) )) 
out <- out %>% mutate(Nas_in_KO_Th17 = apply( out[rownames(meta[meta$condi_real == "IRF4_Ko_1_Th17_D3",])], 1, function(x) sum(is.na(x)) )) 
out <- out %>% mutate(Nas_in_WT_Treg = apply( out[rownames(meta[meta$condi_real == "IRF4_WT_1_TReg_D3",])], 1, function(x) sum(is.na(x)) )) 
out <- out %>% mutate(Nas_in_KO_Treg = apply( out[rownames(meta[meta$condi_real == "IRF4_Ko_1_TReg_D3",])], 1, function(x) sum(is.na(x)) )) 

write.xlsx(out, "output_new/mouse_knockout/MEGA_output.xlsx")

################################################################################
# Heatmaps D3
################################################################################

# palette for heatmaps
# mat_colors <-  colorRampPalette(colors = c(cbPalette[2], cbPalette[4]))(10)
# scales::show_col(mat_colors)
mat_colors <- inferno(20)
# column color annotations!
colors_matrix <- list(
  cell_diff = c("#9C0E0F", "#3E8EB9"),
  IRF4_stat = c("#56B4E9", "#E69F00")
)
names(colors_matrix$cell_diff) <- c("Th17", "Treg")
names(colors_matrix$IRF4_stat) <- c("WT", "KO")
# anno
anno_hm <- meta_complete[c("cell_diff", "IRF4_stat")]
# convert to factor
anno_hm[sapply(anno_hm, is.character)] <- lapply(anno_hm[sapply(anno_hm, is.character)], 
                                       as.factor)
anno_row <- out %>%  dplyr::select("prot.names", "GO_annos_Th17WT_vs_KO", "GO_annos_TregWT_vs_KO")
anno_row <- anno_row %>%  column_to_rownames("prot.names")

experiments_D3 <- experiments_multi[grepl("D3", experiments_multi)]
# set false for filtering NAs instead of imputation
imp = T
data_D3 <- data_complete[c("prot.names", experiments_D3)] %>% column_to_rownames("prot.names")

if (imp != T) {
  # either filter
  data_D3 <- data_D3 %>% filter(!rowSums(is.na(.)) > 0)
} else {
  # or data imputation
  library(impute)
  library(imputeLCMD)
  # imputed <- impute.knn(as.matrix(data_D3), rowmax = 1, colmax=1)
  # data_D3 <- imputed$data %>% as.data.frame()
  data_D3 <- impute.MinDet(as.matrix(data_D3)) %>% as.data.frame()
}



zscores <-  t(scale(t(log2(data_D3))))

hm_inp <- as.matrix(zscores[rownames(data_D3),])
hm <- pheatmap(hm_inp, 
               color=mat_colors,
               cluster_rows=T,
               show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)
ggsave(file = 'output_new/mouse_knockout/heatmap_D3_zscore.pdf',plot = hm,  width=10, height=10)
rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D3_zscore.csv")


# differential genes
signis <- colnames(out)[grepl("signi", colnames(out))]
signis <- signis[2:3] # take only first three
# signis <- signis[4:7] # take last 4
diff_info <- out[c("prot.names", signis)]
# filter significant prots
diff_info <- diff_info %>% filter_at(vars(all_of(signis)), any_vars(. != "ns" & !is.na(.)))


inp <- data_D3 %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% diff_info$prot.names) %>%
  column_to_rownames("prot.names")

zscores <-  t(scale(t(log2(inp))))
hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)
ggsave(file = 'output_new/mouse_knockout/heatmap_D3_differetial_zscore.pdf', plot = hm,  width=10, height=10)
rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D3_differetial_zscore.csv")

# most variable proteins of the differential prots
var_genes <- apply(inp %>%
                     rownames_to_column("prot.names") %>%
                     filter(prot.names %in% diff_info$prot.names) %>%
                     column_to_rownames("prot.names"),
                   1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]

zscores <- inp %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% select_var) %>%
  column_to_rownames("prot.names") %>%
  log2() %>%
  # as.data.frame() %>% 
  mutate(across(everything(), ~ scale(.x)))
hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)

ggsave(file = 'output_new/mouse_knockout/heatmap_D3_diff+mostVar_zscore.pdf', plot = hm,  width=10, height=10)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D3_diff+mostVar_zscore.csv")

# most variable proteins
var_genes <- apply(data_D3,
                   1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:300]

inp <- data_D3 %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% select_var) %>%
  column_to_rownames("prot.names") %>%
  log2() 
zscores <-  t(scale(t(log2(inp))))

hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)

ggsave(file = 'output_new/mouse_knockout/heatmap_D3_mostVar_zscore.pdf', plot = hm,  width=10, height=10)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D3_mostVar_zscore.csv")



logFCs <- str_replace(signis, "signi", "logFC")
diff_log <- c()
diff_log_top100 <- c()
for (i in 1:length(logFCs)) {
  diff_log_top100 <- c(diff_log, out %>%
                  dplyr::arrange(desc(abs(.data[[ logFCs[i] ]]))) %>%
                  head(100) %>%
                  pull(prot.names)
  ) # top 100 each
  diff_log <- c(diff_log, out %>% 
                  filter(abs(.data[[ logFCs[i] ]]) > 1) %>% 
                  pull(prot.names)
  ) # all over logFC1
}

inp <- data_D3 %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% diff_log) %>%
  column_to_rownames("prot.names")

zscores <-  t(scale(t(log2(inp))))
hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)
ggsave(file = 'output_new/mouse_knockout/heatmap_D3_diff+logFC_zscore.pdf',plot = hm,  width=13, height=10)


rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D3_diff+logFC_zscore.csv")

#### add GO + pooling

data_D3_pooled <- data_D3 %>%
  mutate(id = rownames(.)) %>%
  pivot_longer(-id) %>%
  group_by(id, name = str_remove( name, "_[1,2,3]$")) %>%
  summarize(value = mean(value), .groups = "drop") %>% # filter(id == "MPEG1")
  pivot_wider(id_cols = id, names_from = name, values_from = value) %>%
  as.data.frame() %>% 
  column_to_rownames("id")

# anno for pooled data
anno_hm_pool <- anno_hm
anno_hm_pool["new_rn"] <- str_remove(rownames(anno_hm_pool), "_[1-3]$")
anno_hm_pool <- anno_hm_pool %>%  distinct(new_rn, .keep_all = T) 
rownames(anno_hm_pool) <- anno_hm_pool$new_rn
anno_hm_pool <- anno_hm_pool %>% dplyr::select(-new_rn)

# HM

inp <- data_D3_pooled %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% diff_log) %>%
  column_to_rownames("prot.names")

zscores <-  t(scale(t(log2(inp))))

hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm_pool ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)

ggsave(file = 'output_new/mouse_knockout/heatmap_D3_diff+logFC_pooled_zscore.pdf', plot = hm, width=10, height=10)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D3_diff+logFC_pooled_zscore.csv")

# GO Anno
# get clusters and plot with GO overep analysis
# extract clusters
clusts <- sort(cutree(hm$tree_row, k=5))
# get uniprot
df_clusts <- data.frame(prot.names=names(clusts), clusts=clusts)
df_clusts <- df_clusts %>% left_join(out %>% dplyr::select(c("prot.names", "Protein.Group")))

df_clusts["GO_MF"] <- NULL
for (i in 1:length(unique(df_clusts$clusts))){
  prot <- df_clusts %>% filter(clusts == i) %>%  pull(Protein.Group)
  ego <- enrichGO(gene          = prot,
                  universe      = out$Protein.Group,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "UNIPROT",
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  df_clusts[df_clusts$clusts == i, "GO_MF"] <- ego@result$Description[1] 
}
#subset values
df_clusts_ann <- df_clusts %>% column_to_rownames("prot.names") %>% dplyr::select(GO_MF)
# df_clusts_ann$GO_MF <- str_replace_all(df_clusts_ann$GO_MF, ", ", ",\n")
colors_matrix <- list(
  cell_diff = c("#9C0E0F", "#3E8EB9", "#5A5A5A"),
  IRF4_stat = c("#56B4E9", "#E69F00"),
  GO_MF = c("#FD0100", "#F76915", "#EEDE04", "#A0D636", "#2FA236") #, "#333ED4")
)
names(colors_matrix$cell_diff) <- c("Th17", "Treg", "naive")
names(colors_matrix$IRF4_stat) <- c("WT", "KO")
names(colors_matrix$GO_MF) <- unique(df_clusts_ann$GO_MF)

hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm_pool ,
               annotation_colors = colors_matrix,
               annotation_row = df_clusts_ann , 
               annotation_legend = T) 
ggsave(file = 'output_new/mouse_knockout/heatmap_D3_diff+logFC_pooled_zscore+GO.pdf', plot=hm,  width=15, height=10)


rows_txt <- rownames(hm_inp[hm$tree_row$order,])

write.csv(df_clusts[match(rows_txt, df_clusts$prot.names),], "output_new/mouse_knockout/heatmap_D3_diff+logFC_pooled_zscore+GO.csv")

################################################################################
# HMs with d0

colors_matrix <- list(
  cell_diff = c("#9C0E0F", "#3E8EB9", "#5A5A5A"),
  IRF4_stat = c("#56B4E9", "#E69F00")
)
names(colors_matrix$cell_diff) <- c("Th17", "Treg", "naive")
names(colors_matrix$IRF4_stat) <- c("WT", "KO")

data_d0d3 <- data_complete[c("prot.names", experiments_multi)] %>% column_to_rownames("prot.names")

if (imp != T) {
  # either filter
  data_d0d3 <- data_d0d3 %>% filter(!rowSums(is.na(.)) > 0)
} else {
  # or data imputation
  library(impute)
  library(imputeLCMD)
  # imputed <- impute.knn(as.matrix(data_D3), rowmax = 1, colmax=1)
  # data_D3 <- imputed$data %>% as.data.frame()
  data_d0d3 <- impute.MinDet(as.matrix(data_d0d3)) %>% as.data.frame()
}

inp <- data_d0d3 %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% diff_log) %>%
  column_to_rownames("prot.names")

zscores <-  t(scale(t(log2(inp))))

hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)

ggsave(file = 'output_new/mouse_knockout/heatmap_D0D3_diff+logFC_zscore.pdf',plot = hm,  width=13, height=10)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D0D3_diff+logFC_zscore.csv")

# keep only dendro


# build dendro
dend <- hm[[2]] %>%
  as.dendrogram() %>%
  # set("branches_k_color", k=3, colors=c("darkgrey", "darkorange", "steelblue")) %>%
  set("branches_lwd", 2) %>% 
  set("labels_color", "white")

# get annos
lbls <- dend %>% labels
genotype <- ifelse(anno_hm[lbls,]$IRF4_stat == "WT", "#56B4E9", "#E69F00")
cell_type <- ifelse(
  anno_hm[lbls,]$cell_diff == "Th17", "#9C0E0F", ifelse(
    anno_hm[lbls,]$cell_diff == "Treg", "#3E8EB9", "#5A5A5A")
)
clusters <- cutree(hm[[2]], k=3, order_clusters_as_data = FALSE)
bars <- cbind(genotype, cell_type, clusters)
# add anno

pdf(file = 'output_new/mouse_knockout/dendro_D0D3_diff+logFC_zscore.pdf', width=10, height=5)
plot(dend)
colored_bars(colors = bars, y_shift = -0.7)#, dend = dend)
dev.off()

#### pool technical replicates
# cols_norep <- unique(str_remove(colnames(data_d0d3), "_[1,2,3]$"))
data_d0d3_pooled <- data_d0d3 %>%
  mutate(id = rownames(.)) %>%
  pivot_longer(-id) %>%
  group_by(id, name = str_remove( name, "_[1,2,3]$")) %>%
  summarize(value = mean(value), .groups = "drop") %>% # filter(id == "MPEG1")
  pivot_wider(id_cols = id, names_from = name, values_from = value) %>%
  as.data.frame() %>% 
  column_to_rownames("id")



inp <- data_d0d3_pooled %>%
  rownames_to_column("prot.names") %>%
  filter(prot.names %in% diff_log) %>%
  column_to_rownames("prot.names")

zscores <-  t(scale(t(log2(inp))))

hm_inp <- as.matrix(zscores)
hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm_pool ,
               annotation_colors = colors_matrix,
               # annotation_row = anno_row , 
               annotation_legend = T)

ggsave(file = 'output_new/mouse_knockout/heatmap_D0D3_diff+logFC_pooled_zscore.pdf', plot = hm, width=10, height=10)

rows_txt <- rownames(hm_inp[hm$tree_row$order,])
write.csv(rows_txt, "output_new/mouse_knockout/heatmap_D0D3_diff+logFC_pooled_zscore.csv")

# get clusters and plot with GO overep analysis
# extract clusters
clusts <- sort(cutree(hm$tree_row, k=6))
# get uniprot
df_clusts <- data.frame(prot.names=names(clusts), clusts=clusts)
df_clusts <- df_clusts %>% left_join(out %>% dplyr::select(c("prot.names", "Protein.Group")))

df_clusts["GO_MF"] <- NULL
for (i in 1:length(unique(df_clusts$clusts))){
  prot <- df_clusts %>% filter(clusts == i) %>%  pull(Protein.Group)
  ego <- enrichGO(gene          = prot,
                  universe      = out$Protein.Group,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "UNIPROT",
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  df_clusts[df_clusts$clusts == i, "GO_MF"] <- ego@result$Description[1] 
}
#subset values
df_clusts_ann <- df_clusts %>% column_to_rownames("prot.names") %>% dplyr::select(GO_MF)
# df_clusts_ann$GO_MF <- str_replace_all(df_clusts_ann$GO_MF, ", ", ",\n")
colors_matrix <- list(
  cell_diff = c("#9C0E0F", "#3E8EB9", "#5A5A5A"),
  IRF4_stat = c("#56B4E9", "#E69F00"),
  GO_MF = c("#FD0100", "#F76915", "#EEDE04", "#A0D636", "#2FA236", "#333ED4")
)
names(colors_matrix$cell_diff) <- c("Th17", "Treg", "naive")
names(colors_matrix$IRF4_stat) <- c("WT", "KO")
names(colors_matrix$GO_MF) <- unique(df_clusts_ann$GO_MF)

hm <- pheatmap(hm_inp,
               color=mat_colors,
               cluster_rows=T, show_rownames=F, 
               cluster_cols=T,
               annotation_col = anno_hm_pool ,
               annotation_colors = colors_matrix,
               annotation_row = df_clusts_ann , 
               annotation_legend = T) 
ggsave(file = 'output_new/mouse_knockout/heatmap_D0D3_diff+logFC_pooled_zscore+GO.pdf', plot=hm,  width=15, height=10)


rows_txt <- rownames(hm_inp[hm$tree_row$order,])

write.csv(df_clusts[match(rows_txt, df_clusts$prot.names),], "output_new/mouse_knockout/heatmap_D0D3_diff+logFC_pooled_zscore+GO.csv")

