# RNA-Seq Pipeline Documentation

This document provides a detailed overview of a comprehensive RNA-Seq pipeline designed for data processing, alignment, differential expression analysis, and cross-species comparison.

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

The following Docker containers ensure reproducibility and consistency:

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

## Differential Expression Analysis: Retinoblastoma Tumor Study

This study processes RNA-seq data, identifies DEGs, and performs visualization and functional enrichment.

### 1. Data Preparation
Combine **HTSeq** count files:
```r
files <- c("/path/to/file1.htseq.txt", "/path/to/file2.htseq.txt")
combined_df <- lapply(files, read_htseq_file) %>%
  Reduce(function(x, y) merge(x, y, by = c("GeneID", "GeneName")), .)
write.csv(combined_df, "htseq_combined_Counts.csv")
```

### 2. Preprocessing and Filtering
```r
dge <- DGEList(counts = counts, group = factor(c("1", "1", "2", "2", "2")))
dge <- calcNormFactors(dge)
```

### 3. Differential Expression
```r
fit <- glmFit(dge, model.matrix(~group))
lrt <- glmLRT(fit)
DEG_results <- as.data.frame(topTags(lrt, n = Inf))
write.csv(DEG_results, "DEG_results.csv")
```

### 4. Visualization
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

## Cross-Species Analysis: Human ↔ Mouse Orthologs

This section identifies DEGs conserved between human and mouse using ortholog mappings.

### Step 1: Retrieve Orthologs
```r
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
orthologs_human_to_mouse <- getBM(
  attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"),
  filters = "with_mmusculus_homolog",
  values = TRUE, mart = human
)
```

### Step 2: Merge Data
```r
human_to_mouse <- orthologs_human_to_mouse %>%
  inner_join(human_data, by = c("ensembl_gene_id" = "GeneID")) %>%
  inner_join(mouse_data, by = c("mmusculus_homolog_ensembl_gene" = "id"))
```

### Step 3: Visualization
**Venn Diagrams**:
```r
venn.diagram(
  x = list(
    "Human Up" = up_human$ensembl_gene_id,
    "Mouse Up" = up_mouse$ensembl_gene_id
  ), filename = NULL
)
```

**Boxplots**:
```r
ggplot(human_to_mouse, aes(x = direction, y = logFC)) + geom_boxplot()
```

### Step 4: Statistical Tests
```r
wilcox.test(up_human$logFC, human_to_mouse$logFC[human_to_mouse$direction == "Unchanged"])
```

---
