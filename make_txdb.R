# Load required packages
if (!require(AnnotationHub)) BiocManager::install("AnnotationHub")
if (!require(GenomicFeatures)) BiocManager::install("GenomicFeatures")
if (!require(rtracklayer)) BiocManager::install("rtracklayer")
if (!require(biomaRt)) BiocManager::install("biomaRt")

# Set paths
gtf_path <- "hg38/GENCODE/genome/directory/gencode.v44.annotation.gtf"

# Making a TxDb from GTF
txdb <- makeTxDbFromGFF(gtf_path, dataSource = "GENCODE", organism = "Homo sapiens")

# Getting transcript and gene ID columns
keys <- keys(txdb, keytype = "TXNAME")
columns(txdb)
tx2gene <- select(txdb, keys = keys, columns = c("TXID", "TXNAME", "GENEID"), keytype = "TXNAME")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 111)

tx2symbol <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "hgnc_symbol", "gene_biotype"),
                   filters = "ensembl_transcript_id_version",
                   values = tx2gene$TXNAME,
                   mart = ensembl)

# Save references as csvs
write.csv(tx2gene, file = "hg38/GENCODE/genome/directory/GENCODE_tx2gene_v44.csv", row.names = FALSE)
write.csv(tx2symbol, file = "hg38/GENCODE/genome/directory/GENCODE_tx2symbol_v44.csv", row.names = FALSE)
