# HBEC-Media-Comparison

All processing and analysis code for the RNA-seq data associated with the media comparison project, beginning with raw fastq files and concluding with GSEA using DEGs.

### Processing of raw fastq files using Bash

1. Quality of raw files is assessed using fastqc, turned into easily readable html files using multiqc on the fastqc output directory.

2. Remaining adapter sequences and low-quality ends are trimmed from reads using cutadapt.

`cutadapt -q 20
--minimum-length 25
-a AGATCGGAAGAG
-A AGATCGGAAGAG
-o ${ID}_1_adapterAndPhred20Trimmed.fastq.gz
-p ${ID}_2_adapterAndPhred20Trimmed.fastq.gz
${ID}_1_raw.fastq.gz
${ID}_2_raw.fastq.gz`

3. Quality of trimmed files is assessed using fastqc, combined with pre-trimmed output files using multiqc for easy comparison.

4. Trimmed reads are aligned to the transcriptome using STAR, output BAM filed sorted by coordinates.

`STAR --runThreadN 5 --genomeDir hg38/GENCODE/genome/directory --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --outFileNamePrefix ${BASE_FASTQDIR}/${ID}/alignment/${ID}. --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 20 --alignSJoverhangMin  8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin  20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesPrefix ${BASE_FASTQDIR}/${ID}/filtering/ --readFilesIn ${ID}_1_adapterAndPhred20Trimmed.fastq.gz ${ID}_2_adapterAndPhred20Trimmed.fastq.gz`

5. Aligned reads are quantified using Salmon.

`salmon quant --threads 5 --gcBias --gencode -t hg38/GENCODE/genome/directory/gencode.v44.transcripts.fa -l IU -a ${BASE_FASTQDIR}/${ID}/alignment/${ID}.Aligned.toTranscriptome.out.bam -o ${BASE_FASTQDIR}/${ID}/quantification/salmon_quant`

### Differential Expression Analysis using R

1. Make files to map transcript Ensembl ID to gene Ensembl ID and gene symbol (see `make_txdb.R`).

2. Clean data and perform differential expression analysis (see `differential_expression.R`).

3. Using DEG results, perform GSEA (see `GSEA.R`).
