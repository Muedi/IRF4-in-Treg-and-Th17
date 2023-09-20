######################################################################################################
# Script: ChIP-seq Analysis 72h
# Author: Maximilian Sprang, Muedi
# Date: 20.09.2023
# Description: 
# This script overlaps given CHIPseq datasets of Treg and Th17 cells with Enhancer and Silencer Databases
######################################################################################################
# libs
library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)
library(openxlsx)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(stringr)
require(ChIPseeker)
dir.create("output_new/mouse_chip_72h/")

# perm tests are outcommented during writing
perm_test_p <- 110

# load annotation build txdb
txdb <- makeTxDbFromUCSC(genome = "mm10",
                         tablename = "refGene")

# if refGene is not available atm fall back to:
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

anno_mm10 <- toGRanges(txdb, format="GTF", feature="gene")
anno_mm10_biotype <- toGRanges(txdb, format="GTF", feature="gene", colNames="gene_biotype")

# load Enhancer DB 
Enhancers_blood <- read_delim("databases/mm10_lifted.bed", delim = "\t", col_names = FALSE)
colnames(Enhancers_blood) <- c("seqnames", "start", "end", "meta") #, "score", "strand", "cell line", "tissue", "organ", "method
# get cell origins from wacky list
meta_lists <- Enhancers_blood$meta
meta_lists <- str_split(meta_lists, ",")
meta_lists <- lapply(meta_lists, str_extract, ":[A-Z,a-z](.*?)$")
meta_lists <- lapply(meta_lists, str_remove, ":")
meta_lists <- lapply(meta_lists, unique)
unique_origins <- unique(unlist(meta_lists))
Enhancers_blood$meta <- meta_lists
# filtering!!!
# blood related:
# PUER, Natural_killer_cell, Macrophage, CD172+DC, NKT, B_cell_BM, RAW_macrophage
# "CD4+Treg",  "CD8+", CD4+, B_cell_spleen, Th17, "CD4+CD8+", "Neutrophil"
################################################################################
# here we filter for enhancers related to out target cells.                    #
# blood related:                                                               #
# PUER, Natural_killer_cell, Macrophage, CD172+DC, NKT, B_cell_BM, RAW_macrophage
# "CD4+Treg",  "CD8+", CD4+, B_cell_spleen, Th17, "CD4+CD8+", "Neutrophil"     #
# the example below searches for everything related to th17 and Treg!!          #
# Enhancers_blood <- Enhancers_blood %>% filter(grepl("Th17|CD4+Treg", meta) ) #
# the CD4+ example below catches everything that has cd4+ in its description!  #
Enhancers_blood_2 <- Enhancers_blood %>% filter(grepl("CD4+", meta) )          #
################################################################################
# coerce lists to strings
Enhancers_blood_2 <- Enhancers_blood_2 %>% rowwise() %>% 
  mutate(meta = paste(meta, collapse=',')) %>%
  ungroup()
Enhancers_blood_2 <- as.data.frame(Enhancers_blood_2)
Enhancers_blood_2["name"] <- sprintf("Enh_%s",seq(1:dim(Enhancers_blood_2)[1]))
Enhancers_blood_2 <- toGRanges(Enhancers_blood_2)
Enhancers_blood_2 <- unique(Enhancers_blood_2)
# write to xlsx
df_enhancers_blood <- as_tibble(Enhancers_blood_2) %>% mutate(names_enhancers = names(Enhancers_blood_2))
write.xlsx(df_enhancers_blood, "output_new/mouse_chip_72h/enhancers.xlsx")
# load silencer DB
silencers_blood <- read_delim("databases/silencers_Mus_musculus_Blood_lifted_mm10.bed", delim = "\t", col_names = FALSE)
colnames(silencers_blood) <- c("seqnames", "start", "end", "name", "score", "strand", "cell line", "tissue", "organ", "method")
silencers_blood <- as.data.frame(silencers_blood)
rownames(silencers_blood) <- silencers_blood$name
silencers_blood <- toGRanges(silencers_blood)
df_sil_blood <- as_tibble(silencers_blood) %>% mutate(names_enhancers = names(silencers_blood))
write.xlsx(df_sil_blood, "output_new/mouse_chip_72h/silencers.xlsx")

# Overlaps of Th17 vs Treg
data.folder <- "newest_data/chip"
Th17_peaks <- read.xlsx(file.path(data.folder, "USE_EaSeq_v1.2_Th17_windowsize_100_p-val_1E-20_FDR_1E-15_log2_2.1_merge_0.xlsx"))
Th17_peaks["name"] <- sprintf("Th17_%s",seq(1:dim(Th17_peaks)[1]))
Th17_peaks_gr <- toGRanges(Th17_peaks)
Treg_peaks <- read.xlsx(file.path(data.folder, "USE_EaSeq_v1.2_TReg_windowsize_100_p-val_1E-10_FDR_1E-10_log2_1.5_merge_100.xlsx"))
Treg_peaks["name"] <- sprintf("Treg_%s",seq(1:dim(Treg_peaks)[1]))
Treg_peaks_gr <- toGRanges(Treg_peaks)

ol <- findOverlapsOfPeaks(Th17_peaks_gr, Treg_peaks_gr, connectedPeaks = "keepAll")

unique_Treg <- ol$peaklist$Treg_peaks_gr %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup() %>% 
  mutate(peakname = str_extract_all(peakNames, "Treg_\\d+")) %>% 
  unnest(cols=c(peakname))
write.xlsx(unique_Treg, "output_new/mouse_chip_72h/unique_peaks_Treg_bio_vs_rosa.xlsx")

unique_Th17 <- ol$peaklist$Th17_peaks_gr %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup() %>% 
  mutate(peakname = str_extract_all(peakNames, "Th17_\\d+")) %>% 
  unnest(cols=c(peakname))
write.xlsx(unique_Th17, "output_new/mouse_chip_72h/unique_peaks_Th17_bio_vs_rosa.xlsx")


common_peaks <- ol$peaklist$`Th17_peaks_gr///Treg_peaks_gr` %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup() %>% 
  mutate(peakname = str_extract_all(peakNames, "T[a-z, 0-9]*_\\d+")) %>% 
  unnest(cols=c(peakname))
write.xlsx(common_peaks, "output_new/mouse_chip_72h/common_peaks.xlsx")

################################################################################
############################### Th17 ###########################################
################################################################################
##############
############## Enhancer
##############

################################################################################
# Permutation test
# mc.cores can be set to more than 1 in non windows environments only.
# perm_test_p <- peakPermTest(chip_peaks_ranges, Enhancers_blood_2,
#                             TxDb = txdb, mc.cores =1 , ntimes = 1100)
# perm_test_p <- perm_test_p$cntOverlaps$pval
################################################################################
# Venn diagram of overlap.

ol <- findOverlapsOfPeaks(Th17_peaks_gr, 
                          Enhancers_blood_2,
                          connectedPeaks = "keepAll"#, 
                          # maxgap = 100
                          )
pdf("output_new/mouse_chip_72h/Th17_venn_Enhancer_vs_EAseq_peaks_fdr_CD4+.pdf", width = 10, height = 10)
makeVennDiagram(ol, connectedPeaks = "merge",
                NameOfPeaks=c("Th17", "EnhancerDB"),
                fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),
                main="Overlap of lifted CD4+ EnhancerDB and Th17 peaks",
                sub=paste0("p-Value for overlaps: ", perm_test_p)
)
dev.off()


################################################################################
# annotation (not needed)
output_merged <- ol$peaklist$`Th17_peaks_gr///Enhancers_blood_2`

overlaps.anno <- output_merged
ol_df <- as_tibble(overlaps.anno) %>% rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup()

################################################################################
# combine original chip peak file with info of overlaps
# and write to xlsx
ol_df <- ol_df %>% mutate(peakname = str_extract_all(peakNames, "Th17_\\d+"))
ol_df <- ol_df %>% mutate(ols_enhs = str_extract_all(peakNames, "Enh_\\d+")) # overlaps enhancers (sometimes there are multiple)
ol_df <- ol_df %>% unnest(cols=c(peakname))

chip_peaks_combi <- Th17_peaks %>%
  left_join(ol_df %>%
     dplyr::select("peakname",
                   "ols_enhs"), by=c("name"="peakname")) %>% distinct()

##############
############## Silencer
##############

################################################################################
# Permutation test
# mc.cores can be set to more than 1 in non windows environments only.
# perm_test_p <- peakPermTest(chip_peaks_ranges, silencers_blood,
#                             TxDb = txdb, mc.cores=1, ntimes = 1100)
# perm_test_p <- perm_test_p$cntOverlaps$pval

################################################################################
# Venn diagram of overlap.
ol <- findOverlapsOfPeaks(Th17_peaks_gr, silencers_blood, connectedPeaks = "keepAll")
pdf("output_new/mouse_chip_72h/Th17_venn_silencer_vs_EAseq_peaks_fdr.pdf", width = 10, height = 10)
makeVennDiagram(ol,connectedPeaks = "merge",
                NameOfPeaks=c("Th17", "SilencerDB"),
                fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),
                main="Overlap of lifted SilencerDB and Th17 peaks",
                sub=paste0("p-Value for overlaps: ", perm_test_p)
)
dev.off()

output_merged <- ol$peaklist$`Th17_peaks_gr///silencers_blood`

overlaps.anno <- output_merged

ol_df <- as_tibble(overlaps.anno) %>% rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup()

# combine original chip peak file with info of overlaps
ol_df <- ol_df %>% mutate(peakname = str_extract_all(peakNames, "Th17_\\d+"))
ol_df <- ol_df %>% mutate(ols_sils = str_extract_all(peakNames, "MS_\\d+")) # overlaps enhancers (sometimes there are multiple)

ol_df <- ol_df %>% unnest(cols=c(peakname))


chip_peaks_combi <- chip_peaks_combi %>%
                          left_join(ol_df %>% 
                              dplyr::select("peakname",
                                            "ols_sils"), by=c("name"="peakname")) %>% distinct()

##############
############## Triple OL
##############

################################################################################
# Permutation test
# mc.cores can be set to more than 1 in non windows environments only.
# 
# perm_test_p <- peakPermTest(chip_peaks_ranges, silencers_blood,
#                             TxDb = txdb, mc.cores=1, ntimes = 1100)
# perm_test_p <- perm_test_p$cntOverlaps$pval

################################################################################
# Venn diagram of overlap.
ol <- findOverlapsOfPeaks(Th17_peaks_gr, silencers_blood, Enhancers_blood_2, connectedPeaks = "keepAll")
pdf("output_new/mouse_chip_72h/Th17_venn_triple_OL.pdf", width = 10, height = 10)
makeVennDiagram(ol,connectedPeaks = "merge",
                NameOfPeaks=c("Th17", "SilencerDB", "EnhancerDB"),
                fill=c("#CC79A7", "#56B4E9", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2", "#E69F00"), #circle border color
                cat.col=c("#D55E00", "#0072B2", "#E69F00"),
                main="Overlap of lifted SilencerDB, ENhancerDB and Th17 peaks",
                sub=paste0("p-Value for overlaps: ", perm_test_p)
)
dev.off()

output_merged <- ol$peaklist$`Th17_peaks_gr///silencers_blood///Enhancers_blood_2`

overlaps.anno <- output_merged

ol_df <- as_tibble(overlaps.anno) %>% rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup()

# combine original chip peak file with info of overlaps
ol_df <- ol_df %>% mutate(peakname = str_extract_all(peakNames, "Th17_\\d+"))
ol_df <- ol_df %>% mutate(ols_triple = str_extract_all(peakNames, "MS_\\d+|Enh_\\d+")) # overlaps enhancers (sometimes there are multiple)

ol_df <- ol_df %>% unnest(cols=c(peakname))


chip_peaks_combi <- chip_peaks_combi %>%
  left_join(ol_df %>% 
              dplyr::select("peakname",
                            "ols_triple"), by=c("name"="peakname")) %>% distinct()


# annotate
chip_peaks_combi_anno <- ChIPseeker::annotatePeak(GRanges(chip_peaks_combi %>% 
                                                            dplyr::rename("end"="End", "start"="Start")), # quickfix: in Treg end.1 is already renamed...
                                                  TxDb = txdb,
                                                  overlap = "all",
                                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                                  annoDb = "org.Mm.eg.db", 
                                                  tssRegion = c(-3000,3000)) %>% as_tibble()

chip_peaks_combi <- chip_peaks_combi %>%
  left_join(chip_peaks_combi_anno %>%
              dplyr::select("name", 
                            "distanceToTSS",
                            "annotation",
                            "SYMBOL")) %>%
  mutate(unique = ifelse(name %in% unique_Th17$peakname, "unique_Th17", "" ))

# and write to xlsx
write.xlsx(chip_peaks_combi, "output_new/mouse_chip_72h/Th17_ol_EnhancerDB_silencerDB_CD4+.xlsx")


################################################################################
############################### Treg ###########################################
################################################################################

##############
############## Enhancer
##############

# save as bed for uniqpyue peaks
write_delim(Treg_peaks, "output_new/mouse_chip_72h/USE_EaSeq_v1.2_TReg_windowsize_100_p-val_1E-10_FDR_1E-10_log2_1.5_merge_100.bed", delim = "\t", col_names=F)
################################################################################
# Permutation test
# mc.cores can be set to more than 1 in non windows environments only.
# perm_test_p <- peakPermTest(chip_peaks_ranges, Enhancers_blood_2,
#                             TxDb = txdb, mc.cores =1 , ntimes = 1100)
# perm_test_p <- perm_test_p$cntOverlaps$pval

################################################################################
# Venn diagram of overlap.

ol <- findOverlapsOfPeaks(Treg_peaks_gr, Enhancers_blood_2, connectedPeaks = "keepAll")
pdf("output_new/mouse_chip_72h/Treg_venn_Enhancer_vs_EAseq_peaks_fdr_CD4+.pdf", width = 10, height = 10)
makeVennDiagram(ol, connectedPeaks = "merge",
                NameOfPeaks=c("Treg", "EnhancerDB"),
                fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),
                main="Overlap of lifted CD4+ EnhancerDB and Treg peaks",
                sub=paste0("p-Value for overlaps: ", perm_test_p)
)
dev.off()


################################################################################
# annotation (not needed)
output_merged <- ol$peaklist$`Treg_peaks_gr///Enhancers_blood_2`

overlaps.anno <- output_merged

ol_df <- as_tibble(overlaps.anno) %>% rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup()

################################################################################
# combine original chip peak file with info of overlaps
# and write to xlsx
ol_df <- ol_df %>% mutate(peakname = str_extract_all(peakNames, "Treg_\\d+"))
ol_df <- ol_df %>% mutate(ols_enhs = str_extract_all(peakNames, "Enh_\\d+")) # overlaps enhancers (sometimes there are multiple)

ol_df <- ol_df %>% unnest(cols=peakname)

chip_peaks_combi <- Treg_peaks %>%
  left_join(ol_df %>% dplyr::select("peakname",
                                    "ols_enhs"), by=c("name"="peakname")) %>% distinct()

##############
############## Silencer
##############

################################################################################
# Permutation test
# mc.cores can be set to more than 1 in non windows environments only.
# perm_test_p <- peakPermTest(chip_peaks_ranges, silencers_blood,
#                             TxDb = txdb, mc.cores=1, ntimes = 110)
# perm_test_p <- perm_test_p$cntOverlaps$pval

################################################################################
# Venn diagram of overlap.
ol <- findOverlapsOfPeaks(Treg_peaks_gr, silencers_blood, connectedPeaks = "keepAll")
pdf("output_new/mouse_chip_72h/Treg_venn_silencer_vs_EAseq_peaks_fdr.pdf", width = 10, height = 10)
makeVennDiagram(ol,connectedPeaks = "merge",
                NameOfPeaks=c("Treg", "SilencerDB"),
                fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),
                main="Overlap of lifted SilencerDB and Treg peaks",
                sub=paste0("p-Value for overlaps: ", perm_test_p)
)
dev.off()

output_merged <- ol$peaklist$`Treg_peaks_gr///silencers_blood`

overlaps.anno <- output_merged
ol_df <- as_tibble(overlaps.anno) %>% rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup()
# write.xlsx(ol_df, "output_new/mouse_chip_72h/Treg_overlap_silencerDB_fdr.xlsx")

# combine original chip peak file with info of overlaps
ol_df <- ol_df %>% mutate(peakname = str_extract_all(peakNames, "Treg_\\d+"))
ol_df <- ol_df %>% mutate(ols_sils = str_extract_all(peakNames, "MS_\\d+")) # overlaps enhancers (sometimes there are multiple)

ol_df <- ol_df %>% unnest(cols=c(peakname))

chip_peaks_combi <- chip_peaks_combi %>%
  left_join(ol_df %>% 
              dplyr::select("peakname",
                            "ols_sils"), by=c("name"="peakname")) %>% distinct()


##############
############## Triple OL
##############

################################################################################
# Permutation test
# mc.cores can be set to more than 1 in non windows environments only.
# perm_test_p <- peakPermTest(chip_peaks_ranges, silencers_blood,
#                             TxDb = txdb, mc.cores=1, ntimes = 1100)
# perm_test_p <- perm_test_p$cntOverlaps$pval

################################################################################
# Venn diagram of overlap.
ol <- findOverlapsOfPeaks(Treg_peaks_gr, silencers_blood, Enhancers_blood_2, connectedPeaks = "keepAll")
pdf("output_new/mouse_chip_72h/Treg_venn_triple_OL.pdf", width = 10, height = 10)
makeVennDiagram(ol,connectedPeaks = "merge",
                NameOfPeaks=c("Treg", "SilencerDB", "EnhancerDB"),
                fill=c("#CC79A7", "#56B4E9", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2", "#E69F00"), #circle border color
                cat.col=c("#D55E00", "#0072B2", "#E69F00"),
                main="Overlap of lifted SilencerDB, ENhancerDB and Treg peaks",
                sub=paste0("p-Value for overlaps: ", perm_test_p)
)
dev.off()

output_merged <- ol$peaklist$`Treg_peaks_gr///silencers_blood///Enhancers_blood_2`

overlaps.anno <- output_merged

ol_df <- as_tibble(overlaps.anno) %>% rowwise() %>% 
  mutate(peakNames  = paste(peakNames , collapse=',')) %>%
  ungroup()

# combine original chip peak file with info of overlaps
ol_df <- ol_df %>% mutate(peakname = str_extract_all(peakNames, "Treg_\\d+"))
ol_df <- ol_df %>% mutate(ols_triple = str_extract_all(peakNames, "MS_\\d+|Enh_\\d+")) # overlaps enhancers (sometimes there are multiple)

ol_df <- ol_df %>% unnest(cols=c(peakname))


chip_peaks_combi <- chip_peaks_combi %>%
  left_join(ol_df %>% 
              dplyr::select("peakname",
                            "ols_triple"), by=c("name"="peakname")) %>% distinct()



# annotate
chip_peaks_combi_anno <- ChIPseeker::annotatePeak(GRanges(chip_peaks_combi %>% 
                                                            dplyr::rename("end"="End", "start"="Start")), # quickfix
                                                  TxDb = txdb, 
                                                  overlap = "all",
                                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                                  annoDb = "org.Mm.eg.db", 
                                                  tssRegion = c(-3000,3000)) %>% as_tibble()

chip_peaks_combi <- chip_peaks_combi %>%
  left_join(chip_peaks_combi_anno %>%
              dplyr::select("name", 
                            "distanceToTSS",
                            "annotation",
                            "SYMBOL")) %>%
  mutate(unique = ifelse(name %in% unique_Treg$peakname, "unique_Treg", "" ))

# and write to xlsx
write.xlsx(chip_peaks_combi, "output_new/mouse_chip_72h/Treg_ol_EnhancerDB_silencerDB_CD4+.xlsx")


# use gettrack to check for UTR stuff
# library(trackViewer)
# ids <- c("19885")
# symbols <- mget(ids, org.Mm.egSYMBOL)
# geneTrack(ids, txdb, symbols)


# annotation chipseeker
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Mm.eg.db)
library(ChIPpeakAnno)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(openxlsx)
library(GenomicFeatures)
library(tidyverse)

# 
# mm10_gtf <- "C:/Users/masprang/Desktop/Projects/IRF4-project/newest_data/mm10.refGene.gtf.gz"
# txdb <- makeTxDbFromGFF(mm10_gtf)

# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# txdb <- makeTxDbFromUCSC(genome = "mm10",
#                          tablename = "refGene")
###############################################################################
# Th17
################################################################################
data.folder <- "newest_data/chip"
# chip_peaks <- read.xlsx(file.path(data.folder, "USE_EaSeq_v1.2_Th17_windowsize_100_p-val_1E-20_FDR_1E-15_log2_2.1_merge_0.xlsx"))
# chip_peaks_ranges <- toGRanges(chip_peaks)
chip_peaks_ranges <- Th17_peaks_gr
chip_peaks_seeker <- annotatePeak(chip_peaks_ranges,
                                  overlap = "all",
                                  TxDb = txdb, 
                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                  annoDb = "org.Mm.eg.db", 
                                  tssRegion = c(-3000,3000))

write.xlsx(chip_peaks_seeker@anno, "output_new/mouse_chip_72h/Th17_annotation.xlsx")
write.xlsx(chip_peaks_seeker@annoStat, "output_new/mouse_chip_72h/Th17_annotation_stat.xlsx")
# Feature distribution
# upsetplot(chip_peaks_seeker)
pdf(file = 'output_new/mouse_chip_72h/Th17_anno_plot_bar_EAseq-peaks.pdf', width=7, height=3)
plotAnnoBar(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Th17_anno_plot_pie_EAseq-peaks.pdf', width=7, height=7)
plotAnnoPie(chip_peaks_seeker)
dev.off()
# Distribution of TF-binding loci relative to TSS
pdf(file = 'output_new/mouse_chip_72h/Th17_anno_plot_tss_EAseq-peaks.pdf', width=7, height=3)
plotDistToTSS(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Th17_anno_plot_venn_pie_EAseq-peaks.pdf', width=7, height=7)
vennpie(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Th17_anno_upsetplot_EAseq-peaks.pdf', width=7, height=7)
upsetplot(chip_peaks_seeker, vennpie=TRUE)
dev.off()
################################################################################
# Unique peaks
################################################################################

# add information about unique peaks
unique_Th17 <- read.xlsx("output_new/mouse_chip_72h/unique_peaks_Th17_bio_vs_rosa.xlsx")


chip_peaks_ranges <- toGRanges(unique_Th17)
chip_peaks_seeker <- annotatePeak(chip_peaks_ranges,
                                  overlap = "all",
                                  TxDb = txdb, 
                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                  annoDb = "org.Mm.eg.db", 
                                  tssRegion = c(-3000,3000))


write.xlsx(chip_peaks_seeker@anno, "output_new/mouse_chip_72h/Th17_unique_annotation.xlsx")
write.xlsx(chip_peaks_seeker@annoStat, "output_new/mouse_chip_72h/Th17_unique_annotation_stat.xlsx")
# Feature distribution
# upsetplot(chip_peaks_seeker)
pdf(file = 'output_new/mouse_chip_72h/Th17_unique_anno_plot_bar_EAseq-peaks.pdf', width=7, height=3)
plotAnnoBar(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Th17_unique_anno_plot_pie_EAseq-peaks.pdf', width=7, height=7)
plotAnnoPie(chip_peaks_seeker)
dev.off()
# Distribution of TF-binding loci relative to TSS
pdf(file = 'output_new/mouse_chip_72h/Th17_unique_anno_plot_tss_EAseq-peaks.pdf', width=7, height=3)
plotDistToTSS(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Th17_unique_anno_plot_venn_pie_EAseq-peaks.pdf', width=7, height=7)
vennpie(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Th17_unique_anno_upsetplot_EAseq-peaks.pdf', width=7, height=7)
upsetplot(chip_peaks_seeker, vennpie=TRUE)
dev.off()


################################################################################
# Treg
################################################################################
data.folder <- "newest_data/chip"
# chip_peaks <- read.xlsx(file.path(data.folder, "USE_EaSeq_v1.2_TReg_windowsize_100_p-val_1E-10_FDR_1E-10_log2_1.5_merge_100.xlsx"))
# chip_peaks_ranges <- toGRanges(chip_peaks)
chip_peaks_ranges <- Treg_peaks_gr
chip_peaks_seeker <- annotatePeak(chip_peaks_ranges,
                                  overlap = "all",
                                  TxDb = txdb, 
                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                  annoDb = "org.Mm.eg.db", 
                                  tssRegion = c(-3000,3000))

write.xlsx(chip_peaks_seeker@anno, "output_new/mouse_chip_72h/Treg_annotation.xlsx")
write.xlsx(chip_peaks_seeker@annoStat, "output_new/mouse_chip_72h/Treg_annotation_stat.xlsx")
# Feature distribution
# upsetplot(chip_peaks_seeker)
pdf(file = 'output_new/mouse_chip_72h/Treg_anno_plot_bar_EAseq-peaks.pdf', width=7, height=3)
plotAnnoBar(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Treg_anno_plot_pie_EAseq-peaks.pdf', width=7, height=7)
plotAnnoPie(chip_peaks_seeker)
dev.off()
# Distribution of TF-binding loci relative to TSS
pdf(file = 'output_new/mouse_chip_72h/Treg_anno_plot_tss_EAseq-peaks.pdf', width=7, height=3)
plotDistToTSS(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Treg_anno_plot_venn_pie_EAseq-peaks.pdf', width=7, height=7)
vennpie(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Treg_anno_upsetplot_EAseq-peaks.pdf', width=7, height=7)
upsetplot(chip_peaks_seeker, vennpie=TRUE)
dev.off()


################################################################################
# Unique peaks
################################################################################
unique_Treg <-  read.xlsx("output_new/mouse_chip_72h/unique_peaks_Treg_bio_vs_rosa.xlsx")

chip_peaks_ranges <- toGRanges(unique_Treg)
chip_peaks_seeker <- annotatePeak(chip_peaks_ranges,
                                  overlap = "all",
                                  TxDb = txdb, 
                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                  annoDb = "org.Mm.eg.db", 
                                  tssRegion = c(-3000,3000))

write.xlsx(chip_peaks_seeker@anno, "output_new/mouse_chip_72h/Treg_unique_annotation.xlsx")
write.xlsx(chip_peaks_seeker@annoStat, "output_new/mouse_chip_72h/Treg_unique_annotation_stat.xlsx")

# Feature distribution
# upsetplot(chip_peaks_seeker)
pdf(file = 'output_new/mouse_chip_72h/Treg_unique_anno_plot_bar_EAseq-peaks.pdf', width=7, height=3)
plotAnnoBar(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Treg_unique_anno_plot_pie_EAseq-peaks.pdf', width=7, height=7)
plotAnnoPie(chip_peaks_seeker)
dev.off()
# Distribution of TF-binding loci relative to TSS
pdf(file = 'output_new/mouse_chip_72h/Treg_unique_anno_plot_tss_EAseq-peaks.pdf', width=7, height=3)
plotDistToTSS(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Treg_unique_anno_plot_venn_pie_EAseq-peaks.pdf', width=7, height=7)
vennpie(chip_peaks_seeker)
dev.off()
pdf(file = 'output_new/mouse_chip_72h/Treg_unique_anno_upsetplot_EAseq-peaks.pdf', width=7, height=7)
upsetplot(chip_peaks_seeker, vennpie=TRUE)
dev.off()
