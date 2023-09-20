######################################################################################################
# Script: OMICS table
# Author: Maximilian Sprang, Muedi
# Date: 20.09.2023
# Description: 
# This script produces a table containing a condensed version of both the Proteome/Interactome 
# and the ChipSeq datasets.  
######################################################################################################
# libs
library(tidyverse)
library(openxlsx)
library(stringr)
#library(biomaRt)
library(org.Mm.eg.db)


# proteome 
Proteome_Th17vTreg <- read.xlsx("output_new/mouse_72h/results_proteome.xlsx") %>% as_tibble()
remove_x <- colnames(Proteome_Th17vTreg)[grepl("X[0-9][0-9]", colnames(Proteome_Th17vTreg))]
Proteome_Th17vTreg <- Proteome_Th17vTreg %>%
  dplyr::select(-remove_x)

# interactome Wt vs Bio
Interactome_BiovWT <- read.xlsx("output_new/mouse_72h/paired_mega_table_iact.xlsx", startRow = 2) %>%
  as_tibble(.name_repair = "minimal")
names(Interactome_BiovWT) <- make.names(names(Interactome_BiovWT), unique = T)
remove_x <- colnames(Interactome_BiovWT)[grepl("X[0-9][0-9]", colnames(Interactome_BiovWT))]

Interactome_BiovWT <- Interactome_BiovWT %>%
  dplyr::select(-remove_x)

# interactome Cell type vs Cell type
Interactome_Th17vTreg <- read.xlsx("output_new/mouse_72h/paired_results_limma_iact_bioonly.xlsx") %>%
  as_tibble() %>%
  mutate(interactome_Th17vTreg = ifelse(
    logFC > 0.5 & adj.P.Val < 0.01, "th17_enriched", 
      ifelse(
        logFC < -0.5 & adj.P.Val < 0.01, "treg_enriched", "ns")))

# Differentiation experiment
differentiation_exp <- read.xlsx("output_new/mouse_knockout/MEGA_output.xlsx") %>% as_tibble()


# combine tables, info only, no raw intensities
# start with proteome
overview_tab <- Proteome_Th17vTreg %>% 
  dplyr::select(Protein.Group, prot.names, significant_limma) %>%
  dplyr::rename("proteome_Th17vTreg"="significant_limma")
# add interactome bio vs WT
overview_tab <- overview_tab %>%
  full_join(Interactome_BiovWT %>%
              dplyr::select(Protein.Group, prot.names, interactor_Th17, interactor_Treg))

# add interactome cell vs cell
overview_tab <- overview_tab %>%
  full_join(Interactome_Th17vTreg %>%
              dplyr::select(prot.names, Protein.Group, interactome_Th17vTreg )) #, by=c("prot.names", "Protein.Group"))

# add differentiation
signis <- colnames(differentiation_exp)[grepl("signi", colnames(differentiation_exp))]
overview_tab <- overview_tab %>%
  full_join(differentiation_exp %>% dplyr::select(prot.names, Protein.Group, all_of(signis)) )# , by=c("prot.names", "Protein.Group"))
overview_tab <-overview_tab %>% distinct(.keep_all = T)

write.xlsx(overview_tab, "output_new/Overview_omics_no_chip_72h.xlsx")

# map uniprot to ensembl
overview_tab$prot.grp.simple <- word(overview_tab$Protein.Group,1,sep =  ";")

ensembl_ids <- mapIds(org.Mm.eg.db,
                      keys=overview_tab$prot.grp.simple,
                      column = "ENSEMBL",
                      keytype="UNIPROT", 
                      multiVals = 'first')

symbols <- mapIds(org.Mm.eg.db,
                      keys=overview_tab$prot.grp.simple,
                      column = "SYMBOL",
                      keytype="UNIPROT", 
                      multiVals = 'first')

mapped_genes_df <- data.frame(prot.grp.simple=overview_tab$prot.grp.simple,
                              ensembl_id=ensembl_ids,
                              symbol=symbols)

# bind
overview_tab <- overview_tab %>% 
  left_join(mapped_genes_df)

# remove duplicates
overview_tab <-overview_tab %>% distinct(.keep_all = T)
overview_tab_premerge <- overview_tab



# chip
# IRF4 peaks
chip_overlaps_Th17 <- read.xlsx("output_new/mouse_chip_72h/Th17_ol_EnhancerDB_silencerDB_CD4+.xlsx") %>% as_tibble()
chip_overlaps_Treg <- read.xlsx("output_new/mouse_chip_72h/Treg_ol_EnhancerDB_silencerDB_CD4+.xlsx") %>% as_tibble()

overview_tab <- overview_tab %>%
                  left_join(chip_overlaps_Th17 %>%
                             dplyr::select(ols_enhs,
                                           ols_sils,
                                           ols_triple,
                                           distanceToTSS,
                                           annotation,
                                           SYMBOL) %>% 
                             filter(ols_enhs != "" | ols_sils != ""),
                           by=c("symbol"="SYMBOL"),
                           na_matches = "never") %>%
                  dplyr::rename("Th17_peak_ols_enhs"="ols_enhs",
                                "Th17_peak_ols_sils"="ols_sils",
                                "Th17_peak_ols_triple"="ols_triple",
                                "Th17_peak_distanceToTSS"="distanceToTSS",
                                "Th17_peak_annotation"="annotation")

overview_tab <- overview_tab %>%
                  left_join(chip_overlaps_Treg%>%
                              dplyr::select(ols_enhs,
                                            ols_sils,
                                            ols_triple,
                                            distanceToTSS,
                                            annotation,
                                            SYMBOL) %>%
                              filter(ols_enhs != "" | ols_sils != ""),
                            by=c("symbol"="SYMBOL"),
                            na_matches = "never") %>%
                  dplyr::rename("Treg_peak_ols_enhs"="ols_enhs",
                                "Treg_peak_ols_sils"="ols_sils",
                                "Treg_peak_ols_triple"="ols_triple",
                                "Treg_peak_distanceToTSS"="distanceToTSS",
                                "Treg_peak_annotation"="annotation")

# remove duplicate genes and concatenate enhancer/silencer hits
overview_tab <- overview_tab %>%
  group_by(symbol, Protein.Group) %>%
  mutate(Th17_peak_ols_enhs = paste0(Th17_peak_ols_enhs, collapse = ",")) %>%
  mutate(Treg_peak_ols_enhs = paste0(Treg_peak_ols_enhs, collapse = ",")) %>%  
  mutate(Th17_peak_ols_sils = paste0(Th17_peak_ols_sils, collapse = ",")) %>%
  mutate(Treg_peak_ols_sils = paste0(Treg_peak_ols_sils, collapse = ",")) %>%
  mutate(Th17_peak_ols_triple = paste0(Th17_peak_ols_triple, collapse = ",")) %>%
  mutate(Treg_peak_ols_triple = paste0(Treg_peak_ols_triple, collapse = ",")) %>%
  slice_head(n=1) %>% ungroup()
# remove duplicated hits
# this steps removed all Proteins with a NA in their Symbol!
# Could be solved  by grouping with symbol and protname.

# TODO check if still neccessary or need to change


overview_tab <- overview_tab %>%
  mutate(Th17_peak_ols_enhs = strsplit(Th17_peak_ols_enhs, ",") %>%
           map(~toString(unique(.x)))) %>%
  mutate(Th17_peak_ols_sils = strsplit(Th17_peak_ols_sils, ",") %>%
           map(~toString(unique(.x)))) %>%
  mutate(Th17_peak_ols_triple = strsplit(Th17_peak_ols_triple, ",") %>%
           map(~toString(unique(.x)))) %>%
  mutate(Treg_peak_ols_enhs = strsplit(Treg_peak_ols_enhs, ",") %>%
           map(~toString(unique(.x)))) %>%
  mutate(Treg_peak_ols_sils = strsplit(Treg_peak_ols_sils, ",") %>%
           map(~toString(unique(.x)))) %>% 
  mutate(Treg_peak_ols_triple = strsplit(Treg_peak_ols_triple, ",") %>%
           map(~toString(unique(.x)))) 

write.xlsx(overview_tab, "output_new/Overview_omics_72h.xlsx")

