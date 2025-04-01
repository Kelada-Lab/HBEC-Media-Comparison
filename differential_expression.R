# Set seed for reproducibility
set.seed(1)

# Load packages, install if not already
if (!require(tximport)) BiocManager::install("tximport")
if (!require(DESeq2)) BiocManager::install("DESeq2")
if (!require(PCAtools)) BiocManager::install("PCAtools")
if (!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
if (!require(tibble)) install.packages("tibble")
if (!require(tidyverse)) install.packages("tidyverse")

# Set file paths
sampleID_path <- "path/to/metadata.txt"
tx2gene_path <- "hg38/GENCODE/genome/directory/GENCODE_tx2gene_v44.csv"
tx2symbol_path <- "hg38/GENCODE/genome/directory/GENCODE_tx2symbol_v44.csv"
work_path <- "path/to/general/file/directory/"
output_path <- "path/to/output/directory/"

# Read in sample metadata
samples <- read.table(sampleID_path, sep = ",", header = TRUE)
# Scale and center cell type proportions
samples[, c("Ciliated_Prop", "Goblet_Prop")] <- scale(samples[, c("Ciliated_Prop", "Goblet_Prop")], 
                                                      scale = TRUE,
                                                      center = TRUE)

# Import samples and transcript to gene mapping key
files <- file.path("path/to/salmon/output/per/sample/salmon_quant/quant.sf")
names(files) <- paste(samples$Donor, samples$Media, samples$Exposure, samples$Wash_Time, sep = "_")
tx2gene <- read.csv(tx2gene_path, header = TRUE)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Full code should be in the format Media_WashTime_Exposure
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Donor + Full_Code)

# Filter out genes with read counts < 10 in >= 6 samples (one entire category)
keep <- rowSums(counts(ddsTxi) >= 10) >= 6
ddsTxi <- ddsTxi[keep,]

# Run DESeq2
dds <- DESeq(ddsTxi)

# Visualize samples through PCA plot
vst_counts <- assay(vst(dds))
p <- pca(vst_counts, metadata = colData(dds))

biplot(p, 
       showLoadings = FALSE, 
       max.overlaps = Inf,
       colby = "Donor",
       shape = "Media",
       lab = p$metadata$PCA_Label,
       labSize = 5,
       pointSize = 4,
       gridlines.major = FALSE, 
       gridlines.minor = FALSE,
       legendPosition = "none")


# Report results with shrunken fold changes to adjust for genes with low counts in either condition
results_U_2hr_O3vsFA_lfc <- lfcShrink(dds, contrast = c("Full_Code", "U_2hr_O3", "U_2hr_FA"), type="ashr")
results_P_2hr_O3vsFA_lfc <- lfcShrink(dds, contrast = c("Full_Code", "P_2hr_O3", "P_2hr_FA"), type="ashr")

# Convert to dataframe for visualization
df_results_U_2hr_lfc <- as.data.frame(results_U_2hr_O3vsFA_lfc)
df_results_P_2hr_lfc <- as.data.frame(results_P_2hr_O3vsFA_lfc)

# Annotate with gene symbols
tx2symbol <- read.csv(tx2symbol_path, header = TRUE)

df_results_U_2hr_lfc$geneID <- row.names(df_results_U_2hr_lfc)
df_results_U_2hr_lfc <- df_results_U_2hr_lfc %>% left_join(tx2symbol[, c("ensembl_gene_id_version", "hgnc_symbol", "gene_biotype")],
                                                           by = join_by("geneID" == "ensembl_gene_id_version")) %>% 
                                                 relocate(hgnc_symbol, .before = baseMean) %>% 
                                                 distinct()

df_results_P_2hr_lfc$geneID <- row.names(df_results_P_2hr_lfc)
df_results_P_2hr_lfc <- df_results_P_2hr_lfc %>% left_join(tx2symbol[, c("ensembl_gene_id_version", "hgnc_symbol", "gene_biotype")],
                                                           by = join_by("geneID" == "ensembl_gene_id_version")) %>% 
                                                 relocate(hgnc_symbol, .before = baseMean) %>% 
                                                 distinct()

# Annotate DEG direction (or lack of change)
df_results_U_2hr_lfc <- df_results_U_2hr_lfc %>% dplyr::mutate(DE_status = case_when(log2FoldChange >= log2(1.5) & padj <= 0.05 ~ "Up",
                                                                                     log2FoldChange <= -log2(1.5) & padj <= 0.05 ~ "Down",
                                                                                     TRUE ~ "NS"))

df_results_P_2hr_lfc <- df_results_P_2hr_lfc %>% dplyr::mutate(DE_status = case_when(log2FoldChange >= log2(1.5) & padj <= 0.05 ~ "Up",
                                                                                     log2FoldChange <= -log2(1.5) & padj <= 0.05 ~ "Down",
                                                                                     TRUE ~ "NS"))

# Make protein-coding only dataset for main body volcano plots
# Also removing entries that have an Ensembl ID but not a HGNC gene name (likely novel/unconfirmed)
df_results_U_2hr_PC <- df_results_U_2hr_lfc %>% dplyr::filter(gene_biotype == "protein_coding",
                                                              hgnc_symbol != "")

df_results_P_2hr_PC <- df_results_P_2hr_lfc %>% dplyr::filter(gene_biotype == "protein_coding",
                                                              hgnc_symbol != "")

# Volcano Plots
U_2hr_lfc_volcanoPlot_PC <- ggplot(df_results_U_2hr_PC, aes(x = log2FoldChange, y = -log(padj, 10), color = factor(DE_status, levels = c("Down", "NS", "Up")))) +
  geom_point(size = 4, show.legend = FALSE) +
  scale_color_manual(values = c("#619CFF", "gray50", "#F8766D")) +
  geom_vline(xintercept = .58, linetype = "dashed") +
  geom_vline(xintercept = -.58, linetype = "dashed") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dashed") +
  ylim(0, 15.5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"Adj. P")) +
  theme_classic(base_size = 30)

U_2hr_lfc_volcanoPlot_PC

U_2hr_lfc_volcanoPlot_all <- ggplot(df_results_U_2hr_lfc, aes(x = log2FoldChange, y = -log(padj, 10), color = factor(DE_status, levels = c("Down", "NS", "Up")))) +
  geom_point(size = 4, show.legend = FALSE) +
  scale_color_manual(values = c("#619CFF", "gray50", "#F8766D")) +
  geom_vline(xintercept = .58, linetype = "dashed") +
  geom_vline(xintercept = -.58, linetype = "dashed") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dashed") +
  ylim(0, 15.5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"Adj. P")) +
  theme_classic(base_size = 30)

U_2hr_lfc_volcanoPlot_all

P_2hr_lfc_volcanoPlot_PC <- ggplot(df_results_P_2hr_PC, aes(x = log2FoldChange, y = -log(padj, 10), color = factor(DE_status, levels = c("Down", "NS", "Up")))) +
                              geom_point(size = 4, show.legend = FALSE) +
                              scale_color_manual(values = c("gray50")) +
                              geom_vline(xintercept = .58, linetype = "dashed") +
                              geom_vline(xintercept = -.58, linetype = "dashed") +
                              geom_hline(yintercept = -log(0.05, 10), linetype = "dashed") +
                              xlim(-4, 4) +
                              ylim(0, 15.5) +
                              xlab(expression("log"[2]*"FC")) + 
                              ylab(expression("-log"[10]*"Adj. P")) +
                              theme_classic(base_size = 30)

P_2hr_lfc_volcanoPlot_PC

P_2hr_lfc_volcanoPlot_all <- ggplot(df_results_P_2hr_lfc, aes(x = log2FoldChange, y = -log(padj, 10), color = factor(DE_status, levels = c("Down", "NS", "Up")))) +
  geom_point(size = 4, show.legend = FALSE) +
  scale_color_manual(values = c("#619CFF", "gray50", "#F8766D")) +
  geom_vline(xintercept = .58, linetype = "dashed") +
  geom_vline(xintercept = -.58, linetype = "dashed") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dashed") +
  scale_x_continuous(limits = c(-26, 10), breaks = c(-25, -20, -15, -10, -5, 0, 5, 10)) +
  ylim(0, 15.5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"Adj. P")) +
  theme_classic(base_size = 30)

P_2hr_lfc_volcanoPlot_all

# FA exposed media comparison
# Incorporates scaled and centered cell type proportions
# Filter ddsTxi object to only FA samples and change design
ddsTxi_FA <- ddsTxi[,ddsTxi$Exposure == "FA"]
ddsTxi_FA$Media <- factor(ddsTxi_FA$Media, levels = c("U", "P"))
design(ddsTxi_FA) <- ~ Donor + Ciliated_Prop + Goblet_Prop + Media

dds_FA <- DESeq(ddsTxi_FA)

results_FA_2hr_PvsU_lfc <- lfcShrink(dds_FA, contrast = c("Media", "P", "U"), type = "ashr")
df_results_2hr_PvsU_lfc <- as.data.frame(results_FA_2hr_PvsU_lfc)

df_results_2hr_PvsU_lfc$geneID <- row.names(df_results_2hr_PvsU_lfc)
df_results_2hr_PvsU_lfc <- df_results_2hr_PvsU_lfc %>% left_join(tx2symbol[, c("ensembl_gene_id_version", "hgnc_symbol", "gene_biotype")],
                                                                 by = join_by("geneID" == "ensembl_gene_id_version")) %>% 
                                                                 relocate(hgnc_symbol, .before = baseMean) %>% 
                                                                 distinct() %>% 
                                                                 dplyr::mutate(DE_status = case_when(log2FoldChange >= log2(1.5) & padj <= 0.05 ~ "Up",
                                                                                                      log2FoldChange <= -log2(1.5) & padj <= 0.05 ~ "Down",
                                                                                                      TRUE ~ "NS"))

FA_2hr_lfc_volcanoPlot_all <- ggplot(df_results_2hr_PvsU_lfc, aes(x = log2FoldChange, y = -log(padj, 10), color = factor(DE_status, levels = c("Down", "NS", "Up")))) +
                                geom_point(size = 4, show.legend = FALSE) +
                                scale_color_manual(values = c("#00BFC4", "gray50", "#F8766D")) +
                                geom_vline(xintercept = .58, linetype = "dashed") +
                                geom_vline(xintercept = -.58, linetype = "dashed") +
                                geom_hline(yintercept = -log(0.05, 10), linetype = "dashed") +
                                xlim(-10, 10) +
                                xlab(expression("log"[2]*"FC")) + 
                                ylab(expression("-log"[10]*"Adj. P")) +
                                theme_classic(base_size = 30)

FA_2hr_lfc_volcanoPlot_all