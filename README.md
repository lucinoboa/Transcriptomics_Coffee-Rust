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

## Step4: Alignment to the Reference Genome 
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
### Step 5: Quantification of Mapped Reads. 
```bash
featureCounts -p -t exon -g gene_id   -a /data2/lnoboa/ref_genome_coffea/genomic.gtf   -o counts_matrix.txt /data2/lnoboa/roya_transcriptomics/mapping_results/*_sorted.bam
```
