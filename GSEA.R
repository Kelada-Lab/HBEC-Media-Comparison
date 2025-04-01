# Set seed for reproducibility
set.seed(1)

# Load packages, install if not already
if (!require(dplyr)) BiocManager::install("dplyr")
if (!require(RColorBrewer)) BiocManager::install("RColorBrewer")
if (!require(stringr)) BiocManager::install("stringr")
if (!require(clusterProfiler)) BiocManager::install("clusterProfiler")
if (!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
if (!require(janitor)) install.packages("janitor")
if (!require(tidyverse)) install.packages("tidyverse")

### Run enrichGO on UNC DEGs
# Remove genes with NA adj p values (genes that DESeq2 independently filtered out)
df_results_U_2hr_lfc_indFilt <- df_results_U_2hr_lfc[!is.na(df_results_U_2hr_lfc$padj),]

df_results_U_2hr_lfc_indFilt$EnsemblID <- str_split_i(df_results_U_2hr_lfc_indFilt$geneID, "\\.", 1)

# Get DEGs (|log2FC| of 0.58 = |FC| of 1.5)
sig_genes <- df_results_U_2hr_lfc_indFilt[abs(df_results_U_2hr_lfc_indFilt$log2FoldChange) >= 0.58 & df_results_U_2hr_lfc_indFilt$padj <= 0.05,]

# Get background genes and rank by log2FC
geneList_unranked <- df_results_U_2hr_lfc_indFilt[, c("EnsemblID", "log2FoldChange")]
geneList_Ensembl <- geneList_unranked$log2FoldChange
names(geneList_Ensembl) <- geneList_unranked$EnsemblID
geneList_Ensembl <- sort(geneList_Ensembl, decreasing = TRUE)

# Run enrichGO
enrich_results <- enrichGO(gene = sig_genes$hgnc_symbol,
                           universe = names(geneList_Ensembl),
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENSEMBL",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 1,
                           qvalueCutoff = 1,
                           readable = TRUE)

df_enrich_results <- enrich_results@result

# List of prioritized gene sets
selected_genesets <- c("positive regulation of cell death",
                       "response to lipopolysaccharide",
                       "response to molecule of bacterial origin",
                       "cellular response to hypoxia",
                       "cellular response to decreased oxygen levels",
                       "regulation of MAPK cascade",
                       "positive regulation of apoptotic process",
                       "regulation of epithelial cell proliferation",
                       "positive regulation of programmed cell death",
                       "ERK1 and ERK2 cascade",
                       "cellular response to lipid",
                       "response to peptide",
                       "response to oxidative stress",
                       "response to bacterium",
                       "cellular response to lipopolysaccharide",
                       "response to reactive oxygen species",
                       "response to xenobiotic stimulus",
                       "regulation of epithelial cell differentiation",
                       "chemotaxis",
                       "cytokine-mediated signaling pathway",
                       "positive regulation of acute inflammatory response",
                       "prostaglandin secretion",
                       "cellular response to oxidative stress",
                       "regulation of DNA-binding transcription factor activity",
                       "regulation of inflammatory response",
                       "regulation of neutrophil migration",
                       "response to wounding",
                       "cellular response to interleukin-1",
                       "regulation of acute inflammatory response",
                       "cytokine production involved in immune response",
                       "neutrophil chemotaxis")

# Format for discovery and plot making
enrich_results_selected <- df_enrich_results[df_enrich_results$Description %in% selected_genesets,]
enrich_results_selected_short <- enrich_results_selected %>% dplyr::select(Description, geneID) %>% 
  dplyr::mutate(geneID = str_split(geneID, "/"))

# Extract normalized counts and format for plot making
normalized_counts <- DESeq2::counts(dds, normalized = TRUE)

normalized_counts_clean <- normalized_counts %>% as.data.frame() %>% 
  dplyr::mutate(EnsemblID = str_split_i(rownames(.), "\\.", 1)) %>% 
  dplyr::relocate(EnsemblID, .before = Donor1_P_FA_2hr) %>% 
  t() %>% 
  as.data.frame() %>%
  janitor::row_to_names(row_number = 1) %>% 
  dplyr::mutate_all(as.numeric) %>% 
  rownames_to_column(var = "Full_Code") %>% 
  separate(Full_Code, sep = "_", into = c("Donor", "Media", "Exposure", "Wash_Time")) %>% 
  dplyr::mutate_at(c("Donor", "Media", "Exposure"), as.factor)

# Get all genes represented by selected gene sets
all_genes <- unlist(enrich_results_selected_short$geneID)
all_genes <- unique(all_genes)
all_genes <- data.frame("geneID" = all_genes,
                        "numOfSets" = NA,
                        "setNames" = NA)

all_sets <- enrich_results_selected_short$geneID
names(all_sets) <- enrich_results_selected_short$Description

for (i in all_genes$geneID) {
  results <- list.search(all_sets, i %in% .)
  all_genes$numOfSets[all_genes$geneID == i] <- length(results)
  all_genes$setNames[all_genes$geneID == i] <- list(names(results))
}

all_genes <- all_genes[order(all_genes$numOfSets, decreasing = TRUE),]

# Make gene symbol to Ensembl ID reference
symbol_to_EnsemblID <- mapIds(org.Hs.eg.db,
                              keys = all_genes$geneID,
                              column = "ENSEMBL",
                              keytype = "SYMBOL")

# Prep for plotting: build ensembl ID - gene symbol dictionary for selected genes
selected_genes_toPlot <- c("ADAM8", "CXCL2", "CXCL8", "DUSP10", "EDN1", "HMOX1", "IL1B", "MYC", "NR4A1", "ZFP36")

selected_symbolToID <- symbol_to_EnsemblID[names(symbol_to_EnsemblID) %in% selected_genes_toPlot]

selected_counts <- normalized_counts_clean %>% dplyr::select(Donor, Media, Exposure, matches(selected_symbolToID)) %>% 
  dplyr::mutate(Group = paste0(Donor, "_", Media)) %>% 
  relocate(Group, .before = Donor) %>% 
  rename_at(selected_symbolToID, ~names(selected_symbolToID)) %>% 
  pivot_longer(!c(Group, Donor, Media, Exposure), 
               names_to = "GeneName",
               values_to = "NormalizedCount")

plot <- ggplot(selected_counts, aes(x = Exposure, y = NormalizedCount, color = Media, group = Group)) +
  geom_point(aes(shape = Donor), size = 3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~GeneName, nrow = 2, ncol = 5, scales = "free_y") +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  ylab("Normalized Count") +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
plot

### Run gseGO on filtered air PneumaCult vs UNC
df_results_2hr_PvsU_lfc$ensemblID <- str_split_i(df_results_2hr_PvsU_lfc$geneID, "\\.", 1)

for_gsea <- df_results_2hr_PvsU_lfc$log2FoldChange
names(for_gsea) <- df_results_2hr_PvsU_lfc$ensemblID
for_gsea <- sort(for_gsea, decreasing = TRUE)

gsea_results <- gseGO(geneList = for_gsea,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP",
                      minGSSize = 10,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = TRUE)

df_gsea_results <- as.data.frame(gsea_results@result)
PNEU_genesets <- df_gsea_results[df_gsea_results$NES > 0 & df_gsea_results$p.adjust < 0.05,]
UNC_genesets <- df_gsea_results[df_gsea_results$NES < 0 & df_gsea_results$p.adjust < 0.05,]