# RNA-Seq Pipeline
A comprehensive pipeline for RNA-Seq data processing, alignment, differential expression analysis, and cross-species comparison.

## Features
This pipeline supports:
- **Quality Control:** Using FastQC (v0.11.9) and MultiQC (v1.19) for assessing raw read quality.
- **Adapter Trimming:** Using Trim Galore (v0.6.7) to remove adapters and low-quality sequences.
- **Read Alignment:** Using STAR (v2.7.0f) to align reads to reference genomes (GRCh38 or mm10).
- **Quantification:** Using HTSeq (v0.11.2) for gene-level read count generation.
- **Differential Expression Analysis:** Using edgeR (v4.2.2) for identifying differentially expressed genes (DEGs).
- **Gene Ontology (GO) Enrichment:** Using clusterProfiler (v4.12.6) for functional enrichment analysis.
- **Ortholog Mapping and Cross-Species Comparison:** Using biomaRt (v2.60.1) for mapping human-mouse orthologs and comparative analysis.

---

## Docker Containers Used
The following Docker containers ensure a consistent and reproducible environment for each step of the RNA-Seq pipeline:

1. **Base OS**
   - `debian:stable`: Used for running basic scripts.

2. **FastQC**
   - `quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1`: For quality control of sequencing reads.

3. **MultiQC**
   - `staphb/multiqc`: For aggregating FastQC results.
   - `ewels/multiqc:latest`: For generating MultiQC reports.

4. **Trim Galore**
   - `fshiau/trim-galore:latest`: For trimming adapters and low-quality bases.

5. **STAR**
   - `mgibio/star:2.7.0f`: For read alignment and genome indexing.

6. **Samtools**
   - `biocontainers/samtools:v1.9-4-deb_cv1`: For processing and indexing BAM files.

7. **HTSeq**
   - `biocontainers/htseq:v0.11.2-1-deb-py3_cv1`: For gene quantification.

8. **DeepTools**
   - `fshiau/deeptools:latest`: For read coverage analysis.

---

## Commands for Pipeline Steps

### Quality Control
**FastQC:**  
```bash
fastqc -o /path/to/fastqc_results /path/to/trimmed_file.fq.gz
```

**MultiQC:**  
```bash
multiqc /path/to/fastqc_results -o /path/to/multiqc_results
```

---

### Adapter Trimming
**Trim Galore:**  
```bash
trim_galore -j 8 --gzip -o "${output_dir}" --paired --clip_R2 4 "$r1" "$r2"
```

---

### Genome Indexing
**STAR (Generate Genome Index):**  
```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ${genome_dir} \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v43.primary_assembly.annotation.gtf \
     > star_genome_generate.log 2>&1
```

---

### Read Alignment
**STAR (Align Reads):**  
```bash
STAR --genomeDir "${genome_dir}" \
     --runThreadN 16 \
     --genomeLoad LoadAndKeep \
     --readFilesCommand zcat \
     --readFilesIn "${r1}" "${r2}" \
     --outFileNamePrefix "${output_dir}/${sample_base}_" \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 20000000000
```

---

### BAM Processing
**Filter, Sort, and Index with Samtools:**  
```bash
samtools view -F 1548 -q 30 -bSu "$bam" | samtools sort -o "$cleansort"
```

---

### Gene Quantification
**HTSeq-count:**  
```bash
htseq-count -r pos -s no --additional-attr=gene_name -f bam $cleansort \
/path/to/gencode.v43.primary_assembly.annotation.gtf > $outcount
```

---

### Coverage Analysis
**bamCoverage:**  
```bash
bamCoverage -b "$outbamblk" -o "$coverage" \
    -p max \
    --binSize 10 \
    -e \
    -bl /path/to/.bed \
    --normalizeUsing CPM
```

---

# Human Tumor Analysis: Retinoblastoma

This analysis involves RNA-seq data processing, differential expression analysis, visualization, and pathway enrichment for human retinoblastoma tumors.

---

## Key Steps in Analysis

### 1. **Data Preparation**
#### Combine HTSeq Count Files
```r
# List and read HTSeq files
files <- c("/path/to/file1.htseq.txt", "/path/to/file2.htseq.txt", ...)
combined_df <- lapply(files, read_htseq_file) %>%
  Reduce(function(x, y) merge(x, y, by = c("GeneID", "GeneName")), .)

# Save combined count data
write.csv(combined_df, "htseq_combined_Counts.csv", row.names = FALSE)
```

---

### 2. **Preprocessing and Filtering**
#### Load Count Data and Filter
```r
# Load count data
counts <- read.csv("htseq_combined_Counts.csv", header = TRUE)
counts <- counts[!duplicated(counts[, 2]), ]  # Remove duplicates
row.names(counts) <- counts[, 2]
counts <- counts[, -c(1, 2)]

# Create DGEList and filter low-count genes
dge <- DGEList(counts = counts, group = factor(c("1", "1", "2", "2", "2")))
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
```

---

### 3. **Differential Expression Analysis**
#### Perform Analysis
```r
# Estimate dispersions and fit GLM model
dge <- estimateDisp(dge)
fit <- glmFit(dge, model.matrix(~group))
lrt <- glmLRT(fit)

# Extract DE genes
DEG_results <- as.data.frame(topTags(lrt, n = Inf))
write.csv(DEG_results, "DEG_results_genename_fc.csv", row.names = TRUE)
```

---

### 4. **Visualization**
#### Volcano Plot
```r
library(ggplot2)
volcano_data <- data.frame(Gene = rownames(DEG_results), logFC = DEG_results$logFC, negLogP = -log10(DEG_results$PValue))
volcano_data$Status <- "Not Significant"
volcano_data$Status[volcano_data$logFC > 1 & volcano_data$FDR <= 0.01] <- "Upregulated"
volcano_data$Status[volcano_data$logFC < -1 & volcano_data$FDR <= 0.01] <- "Downregulated"

ggplot(volcano_data, aes(x = logFC, y = negLogP, color = Status)) +
  geom_point() + scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-Value")
```

#### Heatmap
```r
library(pheatmap)
log_cpm <- cpm(dge, log = TRUE)
z_scores <- t(scale(t(log_cpm[rownames(dge) %in% top_genes, ])))

pheatmap(z_scores, main = "Top DE Genes Heatmap", cluster_rows = TRUE, cluster_cols = TRUE)
```

---

### 5. **Pathway Enrichment**
#### GO and KEGG Enrichment
```r
library(clusterProfiler)
significant_genes <- rownames(DEG_results)[abs(DEG_results$logFC) >= 1 & DEG_results$FDR <= 0.01]
entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG
kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.1)
dotplot(kegg, showCategory = 20)

# GO
go <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.1)
dotplot(go, showCategory = 20)
```

---
