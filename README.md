# Transcriptomics_Coffee-Rust
This is a brief guide on how to carry out the quality control, trimming, alignment, and differential expression analysis of RNAseq of coffee plants inoculated with the pathogen Hemileia vastatrix. 

## Step 1: Quality Control of Raw Data 
### Navigate to the folder that contains the raw ".fq.gz" files.  
```bash
cd /data2/lnoboa/roya_transcriptomics/data_cafe_roya/raw_data_roya
```
### Activate the conda environment. 
```bash
conda activate ______________
```
### Run FastQC to assess the quality of the sequences.  
```bash
fastqc *.fq.gz -o /data2/lnoboa/roya_transcriptomics/fastqc_cafe_roya/fastqc_rawdata
```
### Run MultiQC in the same folder. 
```bash
multiqc .
```
[Check the report](multiqc_report_roya_raw_data.html)

## Step 2: Trimming 
###Run the Trimmomatic tool to remove low quality reads and adapters. 
```bash
TRIMMOMATIC=/home/jupyter-alumno7/.conda/envs/rnseq/bin/trimmomatic
ADAPTERS=/home/jupyter-alumno7/.conda/envs/rnseq/share/trimmomatic/adapters/TruSeq3-PE.fa
THREADS=4
INPUT_DIR=/data2/lnoboa/tricho_transcriptomics/raw_data_tricho
OUTPUT_DIR=/data2/lnoboa/tricho_transcriptomics/trimmed_tricho

for f1 in $INPUT_DIR/*_1.fq.gz; do
    base=$(basename "$f1" _1.fq.gz)
    f2="$INPUT_DIR/${base}_2.fq.gz"
     
    $TRIMMOMATIC PE -threads $THREADS \
        "$f1" "$f2" \
        "$OUTPUT_DIR/${base}_1.trim.fq.gz" "$OUTPUT_DIR/${base}_1un.trim.fq.gz" \
        "$OUTPUT_DIR/${base}_2.trim.fq.gz" "$OUTPUT_DIR/${base}_2un.trim.fq.gz" \
        ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
  done
```

## Step 3: Quality Control of the Trimmed Sequences 
### Navigate to the folder that contains the trimmed ".fq.gz" files.  
```bash
cd /data2/lnoboa/roya_transcriptomics/data_cafe_roya/data_trimmed
```
### Activate the conda environment. 
```bash
conda activate ______________
```
### Run FastQC to assess the quality of the sequences.  
```bash
fastqc *.fq.gz -o /data2/lnoboa/roya_transcriptomics/fastqc_cafe_roya/fastqc_trimmed
```
### Run MultiQC in the same folder. 
```bash
multiqc .
```
[Check the report](multiqc_report_roya_trimmed.html)

## Step 4: Alignment to the Reference Genome 
### Create a folder for results. 
```bash
cd data2/lnoboa/roya_transcriptomics/
mkdir mapping_results
```
### Upload the genome index. 
[Reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036785885.1/)
```bash
cd /data2/lnoboa/ref_genome_coffea
hisat2-build -p 4 GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna coffea_index
```
### Mapping of sequences. 
#### Alignment of the trimmed sequences to the reference genome. 
```bash
cd /data2/lnoboa/roya_transcriptomics/data_cafe_roya/data_trimmed
for sample in *_1.trim.fq.gz; do
    base=$(basename $sample _1.trim.fq.gz)
    hisat2 -p 4 \
        -x /data2/lnoboa/ref_genome_coffea/coffea_index \
        -1 ${base}_1.trim.fq.gz \
        -2 ${base}_2.trim.fq.gz \
        -S /data2/lnoboa/mapping_results/${base}.sam
done
```
### Convert SAM to BAM. 
```bash
for f in *.sam; do
    base=$(basename "$f" .sam)
    samtools sort "$f" -o "${base}_sorted.bam"
```
## Step 5: Quantification of Mapped Reads. 
```bash
featureCounts -p -t exon -g gene_id   -a /data2/lnoboa/ref_genome_coffea/genomic.gtf   -o counts_matrix.txt /data2/lnoboa/roya_transcriptomics/mapping_results/*_sorted.bam
```

## Step 6: Differential Expression Analysis with DESeq2
### Install and load the recquired packages. 
```bash
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
install.packages(c("ggplot2", "pheatmap"))
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
```
### Set working directory. 
```r
setwd("D:/lucianoboa/royatranscriptomics/analysis/featureCounts")
```
### Load the count matrix. 
```r
countData <- read.table("counts_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
```
### Check column names.
```r
colnames(countData)
```
### Select only the count columns (columns 6 to 19).
```r
countData <- countData[, 6:19]
```
### Check dimensions of the filtered count matrix.
```r
dim(countData)
```
```r
[1] 75881    14
```
### Rename columns with simpler sample names.
```r
colnames(countData) <- c("H10","H11","H12","H13","H14","H15","H16","H9",
                         "T1","T4","T5","T6","T7","T8")
```
### Define experimental conditions.
#### Eight samples H = 24 hours, Six samples T = 0 hours.
```r
condition <- factor(c(rep("H24", 8), rep("T0", 6)))
```
### Create colData with experimental metadata.
```r
colData <- data.frame(row.names = colnames(countData),
                      condition = condition)
```
### Assembly of the differential expression matrix.
```r
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "H24", "T0"))
res <- res[order(res$padj), ]
head(res)
```
### Top differentially expressed genes (DEGs)
```r
head(res)
```
#### The following shows the top six genes sorted by adjusted p-value from the DESeq2 analysis (H24 vs T0):
```r
log2 fold change (MLE): condition H24 vs T0 
Wald test p-value: condition H24 vs T0 
DataFrame with 6 rows and 6 columns
              baseMean log2FoldChange     lfcSE      stat      pvalue        padj
LOC113711922  1743.061        5.70644  0.309049   18.4645 3.98535e-76 1.62539e-71
LOC113695446  2973.581        6.29184  0.370363   16.9883 1.00210e-64 1.44129e-60
LOC113706933  1380.963        4.13309  0.243337   16.9850 1.06019e-64 1.44129e-60
LOC113735746  1816.273       -1.73424  0.112012  -15.4827 4.54336e-54 4.63241e-50
LOC113737160   636.291        2.50899  0.170068   14.7529 2.94895e-49 2.40540e-45
LOC113695056  8194.214        5.60521  0.382360   14.6595 1.17121e-48 7.96114e-45
```
### DESeq2 Summary of Differential Expression (H24 vs T0)
```r
summary(res)
```
#### The summary below shows the overall statistics from the DESeq2 analysis:
```r
Summary of DESeq2 results:
Total genes with non-zero counts: 47,227

Genes with adjusted p-value < 0.1:
  - Upregulated (LFC > 0): 6,277 (13%)
  - Downregulated (LFC < 0): 6,023 (13%)

Other filtering notes:
  - Outliers [1]: 302 (0.64%)
      [1] See 'cooksCutoff' argument in ?results
  - Low counts [2]: 6,438 (14%) (mean count < 2)
      [2] See 'independentFiltering' argument in ?results
```

[PCA Plot](plotPCA_roya-transcriptomics.png)
