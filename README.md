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

