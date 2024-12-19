# RNA-Seq Pipeline Documentation

This project uses an RNA-Seq pipeline to:

**Study Human Retinoblastoma Tumors:** Analyze gene expression differences between tumor and normal tissues.

**Analyze Cross-Species Gene Expression:** Compare gene expression between human and mouse to explore species-specific patterns.

---

## Pipeline Features

The pipeline incorporates the following steps:

- **Quality Control**: Utilizing **FastQC** and **MultiQC** for assessing raw read quality.
- **Adapter Trimming**: Removing adapters and low-quality sequences using **Trim Galore**.
- **Read Alignment**: Aligning reads to reference genomes (**GRCh38** or **mm10**) using **STAR**.
- **Gene Quantification**: Generating gene-level counts with **HTSeq**.
- **Differential Expression Analysis**: Identifying DEGs using **edgeR**.
- **Functional Enrichment**: Performing **Gene Ontology (GO)** and **KEGG** enrichment using **clusterProfiler**.
- **Cross-Species Comparison**: Mapping orthologs between human and mouse with **biomaRt**.

---

## Docker Containers Utilized

Docker containers:

| Step                 | Docker Container                                   |
|----------------------|---------------------------------------------------|
| Base OS             | `debian:stable`                                   |
| Quality Control     | `quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1` |
| MultiQC             | `staphb/multiqc`, `ewels/multiqc:latest`          |
| Adapter Trimming    | `fshiau/trim-galore:latest`                       |
| Read Alignment      | `mgibio/star:2.7.0f`                              |
| BAM Processing      | `biocontainers/samtools:v1.9-4-deb_cv1`           |
| Gene Quantification | `biocontainers/htseq:v0.11.2-1-deb-py3_cv1`       |
| Coverage Analysis   | `fshiau/deeptools:latest`                         |

---

## Workflow Commands

### Step 1: Quality Control
**FastQC**:
```bash
fastqc -o /path/to/fastqc_results /path/to/trimmed_file.fq.gz
```
**MultiQC**:
```bash
multiqc /path/to/fastqc_results -o /path/to/multiqc_results
```

---

### Step 2: Adapter Trimming
**Trim Galore**:
```bash
trim_galore -j 8 --gzip -o "${output_dir}" --paired "$r1" "$r2"
```

---

### Step 3: Read Alignment
**STAR Genome Indexing**:
```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ${genome_dir} \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v43.primary_assembly.annotation.gtf
```
**STAR Alignment**:
```bash
STAR --genomeDir "${genome_dir}" \
     --runThreadN 16 \
     --readFilesCommand zcat \
     --readFilesIn "${r1}" "${r2}" \
     --outSAMtype BAM SortedByCoordinate
```

---

### Step 4: Gene Quantification
**HTSeq**:
```bash
htseq-count -r pos -s no --additional-attr=gene_name -f bam $cleansort \
/path/to/gencode.v43.primary_assembly.annotation.gtf > $outcount
```

---

### Step 5: Coverage Analysis
**bamCoverage**:
```bash
bamCoverage -b "$outbamblk" -o "$coverage" \
    -p max --binSize 10 --normalizeUsing CPM
```

---

---

## Differential Expression Analysis: Retinoblastoma Tumor Study

### **Which genes are differentially expressed in retinoblastoma tumors compared to normal tissues?**

This study processes RNA-seq data to identify differentially expressed genes (DEGs) in retinoblastoma tumors and employs visualization techniques to present the findings.

### 1. Data Preparation
The **HTSeq** count files are combined for downstream analysis:
```r
files <- c("/path/to/file1.htseq.txt", "/path/to/file2.htseq.txt")
combined_df <- lapply(files, read_htseq_file) %>%
  Reduce(function(x, y) merge(x, y, by = c("GeneID", "GeneName")), .)
write.csv(combined_df, "htseq_combined_Counts.csv")
```

### 2. Preprocessing and Filtering
Normalization and filtering ensure data quality:
```r
dge <- DGEList(counts = counts, group = factor(c("1", "1", "2", "2", "2")))
dge <- calcNormFactors(dge)
```

### 3. Differential Expression
DEGs are identified through statistical modeling:
```r
fit <- glmFit(dge, model.matrix(~group))
lrt <- glmLRT(fit)
DEG_results <- as.data.frame(topTags(lrt, n = Inf))
write.csv(DEG_results, "DEG_results.csv")
```

### 4. Visualization
Key visualizations to reveal patterns in the data:

**Volcano Plot**:
```r
ggplot(volcano_data, aes(x = logFC, y = negLogP, color = Status)) +
  geom_point() + labs(title = "Volcano Plot")
```

**Heatmap**:
```r
pheatmap(log_cpm, main = "Top DE Genes Heatmap")
```

---

### **What are the functional roles of the identified differentially expressed genes, and how do they contribute to specific pathways?**

Pathway enrichment analyses provide insights into the biological significance of DEGs.

**GO and KEGG Enrichment**:
```r
library(clusterProfiler)
significant_genes <- rownames(DEG_results)[abs(DEG_results$logFC) >= 1 & DEG_results$PValue <= 0.1]
entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG Enrichment
kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.1)
dotplot(kegg, showCategory = 20)

# GO Enrichment
go <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.1)
dotplot(go, showCategory = 20)
```

---
## Rod and Cone Gene Expression Analysis

### **What are the overall expression trends for rod- and cone-enriched genes across normal (N) and tumor (T) samples, as observed through RNA-Seq data?**
---

#### **Overview of the Approach**

1. **Gene Selection**: Predefined lists of rod- and cone-enriched genes were used to filter the RNA-Seq dataset.
2. **Normalization**: Log-transformed counts per million (CPM) and Z-score normalization were applied to standardize expression data across samples.
3. **Visualization**: Heatmaps and boxplots were generated to compare expression patterns and summarize trends between normal and tumor samples.

---

#### **Key Steps and Code Snippets**

**1. Filtering Rod- and Cone-Enriched Genes**

Rod- and cone-enriched genes were selected from the dataset using predefined gene lists:

```r
# Define rod- and cone-enriched gene lists
rod_genes <- toupper(c("SAMD11", "RHO", "PDE6A", "GNAT1", "GNGT1", "CNGB1"........))
cone_genes <- toupper(c("GNAT2", "CNGB3", "PDE6H", "GUCA1C", "OPN1MW"........))

# Filter counts for rod and cone genes
counts$GeneName <- toupper(counts$GeneName)  # Ensure consistency
rod_gene_counts <- counts[counts$GeneName %in% rod_genes, ]
cone_gene_counts <- counts[counts$GeneName %in% cone_genes, ]
```

---

**2. Normalizing Data and Generating Heatmaps**

Log-transformed CPM values were computed and standardized using Z-scores, followed by heatmap visualization:

```r
# Compute log-transformed CPM and Z-scores
log_cpm <- cpm(dge, log = TRUE)
z_scores <- t(scale(t(log_cpm)))

# Generate heatmap for rod genes
pheatmap(z_scores, 
         main = "Rod Gene Expression Heatmap", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE, 
         cluster_cols = TRUE)
```

---

**3. Comparing Average Gene Expression**

Boxplots were created to summarize average expression differences between normal and tumor samples:

```r
# Calculate average expression for rods
rod_gene_counts$N_avg <- rowMeans(rod_gene_counts[, c("N1", "N2")])
rod_gene_counts$T_avg <- rowMeans(rod_gene_counts[, c("T1", "T2", "T3")])

# Melt data for visualization
rod_melt <- melt(rod_gene_counts[, c("N_avg", "T_avg")], 
                 variable.name = "Condition", 
                 value.name = "Expression")

# Plot boxplot
ggplot(rod_melt, aes(x = Condition, y = log2(Expression + 1), fill = Condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(title = "Rod Gene Expression Between Conditions",
       x = "Condition", y = "log2(Average Expression + 1)")
```

---
## Cross-Species Analysis: Human ↔ Mouse Orthologs

## **Are the observed gene expression changes in retinoblastoma conserved between human and mouse models?***

This analysis aims to identify and visualize conserved gene expression changes in retinoblastoma between human and mouse models using ortholog mappings and differential expression analysis.

---

#### **Step 1: Retrieve Orthologs**

Retrieve orthologs between human and mouse genomes using the Ensembl database:

```r
library(biomaRt)
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
orthologs <- getBM(
  attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"),
  filters = "with_mmusculus_homolog",
  values = TRUE,
  mart = human,
  uniqueRows = TRUE
)
```

### **Step 2: Merge Data**

Combine human differential expression data with orthologs and subsequently merge with mouse expression data to create a unified dataset:

```r
# Preprocess and clean Ensembl IDs
human_data$GeneID <- gsub("\\..*", "", human_data$GeneID)
orthologs$ensembl_gene_id <- gsub("\\..*", "", orthologs$ensembl_gene_id)
mouse_data$id <- gsub("\\..*", "", mouse_data$id)

# Merge human data with orthologs
merged_human_orthologs <- merge(human_data, orthologs, by.x = "GeneID", by.y = "ensembl_gene_id")

# Merge with mouse data
merged_human_mouse <- merge(merged_human_orthologs, mouse_data, 
                            by.x = "mmusculus_homolog_ensembl_gene", by.y = "id")
```

#### **Step 3: Classify Genes by Expression**

Classify genes based on fold-change and significance thresholds:

```r
upregulated_human <- merged_human_mouse$ensembl_gene_id[merged_human_mouse$logFC >= 1 & merged_human_mouse$FDR < 0.01]
downregulated_human <- merged_human_mouse$ensembl_gene_id[merged_human_mouse$logFC <= -1 & merged_human_mouse$FDR < 0.01]

upregulated_mouse <- merged_human_mouse$ensembl_gene_id[merged_human_mouse$tcKO_v_CreNeg_LogFC >= 1 & merged_human_mouse$tcKO_v_CreNeg_FDR < 0.01]
downregulated_mouse <- merged_human_mouse$ensembl_gene_id[merged_human_mouse$tcKO_v_CreNeg_LogFC <= -1 & merged_human_mouse$tcKO_v_CreNeg_FDR < 0.01]
```

#### **Step 4: Visualize Results**

- **Venn Diagrams**:
  
Visualize overlaps between up- and down-regulated orthologs in human and mouse:

```r
venn.diagram(
  x = list(
    "Human Upregulated" = upregulated_human,
    "Mouse Upregulated" = upregulated_mouse
  ),
  filename = NULL,
  main = "Overlap of Upregulated Genes"
)
```

- **Boxplots**:
  
Compare fold changes between human and mouse orthologs:

```r
ggplot(merged_human_mouse, aes(x = direction, y = tcKO_v_CreNeg_LogFC, fill = direction)) +
  geom_boxplot(outlier.color = "yellow", notch = TRUE, width = 0.5, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(title = "Gene Expression Across Human and Mouse Orthologs",
       x = "Expression Category (Human)", y = "Mouse Fold Change")
```

#### **Step 5: Statistical Tests**

Perform Wilcoxon tests to assess differences in fold changes:

```r
wilcox.test
```

---

