
<img width="955" alt="Screenshot 2024-09-09 at 09 02 30" src="https://github.com/user-attachments/assets/12edca2b-af13-4d28-ab37-133a8fe49ea4">

This workflow combines RNA-seq data analysis, differential gene expression identification, enrichment analysis, and miRNA-target gene prediction to provide comprehensive insights into the molecular mechanisms involved in the study.

- RNA-seq Data Processing: Raw RNA-seq data undergoes quality control (FastQC) and read quantification (featureCounts) to measure gene expression levels.

- Differential Gene Expression (DEG) Analysis: The identified gene expression data is used to detect differentially expressed genes (DEGs) between different experimental conditions using statistical analysis.

- miRNA and Gene Target Prediction: After DEG analysis, miRNA-target prediction is conducted to identify miRNAs that may regulate the expression of the DEGs. This helps in understanding post-transcriptional regulatory mechanisms.

- Gene Ontology (GO) Enrichment Analysis: Enrichment analysis is performed on the DEGs to identify over-represented biological processes, molecular functions, and cellular components, providing context for the functions of the genes identified.

- Topology Analysis: Topological analysis is performed on the network of DEGs (both up-regulated, down-regulated, and all). This step is crucial for identifying key genes (hubs) in the biological network and analyzing the connectivity and influence of these genes in the overall regulatory network. Hubs in these networks are often biologically significant as they may represent genes critical for specific cellular functions or disease mechanisms. Bottleneck nodes may also be identified, which play important roles in connecting different modules or pathways in the network.





## fqs_to_featurecounts:
```
python3 rnaseq_pipeline.py --fastq-dir /path/to/fastqs \
                           --index-dir /path/to/star_index \
                           --gtf-file /path/to/annotation.gtf \
                           --output-dir /path/to/output \
                           --threads 8
```

Assumptions:
- Paired-end FASTQs are named with _R1 and _R2 in their filenames.
- Single-end FASTQs only contain _R1.
- All FASTQs for different samples are located in the same directory (--fastq-dir).
## or
```
snakemake -s rnaseq_pipeline.smk --cores 8
```
- SAMPLES, FASTQ1, FASTQ2: Extracts the sample names and FASTQ paths using glob_wildcards.
- rule star_align: Runs STAR to align paired-end FASTQ files to the reference genome and produces sorted BAM files and read counts.
- rule featurecounts: Runs featureCounts to count gene-level features from the aligned BAM files.
- rule all: This is the top-level rule that specifies the final outputs (all BAM files and gene count files) so Snakemake knows the dependencies.


Folder Structure

```
project_directory/
├── data/
│   ├── samplesheet.csv  # Sample sheet with FASTQ file paths
│   ├── fastq/  # Directory containing FASTQ files
├── reference/
│   ├── genome_index/  # STAR index files
│   ├── annotation.gtf  # GTF annotation file
├── results/
│   ├── alignments/  # STAR output
│   ├── counts/  # featureCounts output
├── rnaseq_pipeline.smk  # Snakemake workflow
└── config.yaml  # Configuration file
```

## deg: 
### deg between two conditions:
``` 
Rscript deg_analysis.R --raw.reads.csv "featureCount_two_cond.csv" \
                       --colname.case "case_a1" \
                       --number.case.samples 4 \
                       --named.case "case_a" \
                       --colname.control "case_b1" \
                       --number.control.samples 4 \
                       --named.control "case_b" \
                       --number.logFC 1.5 \
                       --number.FDR 0.01
```

Command-Line Arguments:
- raw.reads.csv: The CSV file containing raw counts.
- colname.case: Column name for the first case sample.
- number.case.samples: Number of case samples.
- named.case: Name for the case group.
- colname.control: Column name for the first control sample.
- number.control.samples: Number of control samples.
- named.control: Name for the control group.
- number.logFC: The log fold change threshold.
- number.FDR: The False Discovery Rate threshold.

### deg between two paried-sample conditions:
```
Rscript deg_pd_analysis.R --count-data-file "featureCount_pd_cond.csv" \
                          --output-prefix "vehicle_vs_evhd" \
                          --group-labels "Group1,Group2,Group3,Group4" \
                          --condition-labels "c1,c2" \
                          --n-replicates-per-group 2 \
                          --n-groups 4 \
                          --ref-condition "c1" \
                          --target-condition "c2"  
```
- count-data-file "data_pTreg.csv": Path to CSV with gene counts (rows: genes, columns: samples).
- output-prefix "c1_vs_c2": Prefix for output files.
- group-labels "Group1,Group2,Group3,Group4": Comma-separated list of group labels.
- condition-labels "c1,c2": Comma-separated list of condition labels.
- n-replicates-per-group 2: Number of replicates per group.
- n-groups 4: Number of groups in the dataset.
- ref-condition "c1": Reference condition for comparison.
- target-condition "c2": Condition to compare against the reference.

## perform enrichment analysis
```
Rscript perform_enrichment_analysis.R file_name 3702
```
- file_name: The exact name of the CSV file you want to process e.g. perform_enrichment.csv.
- organism_id: The organism ID used for the enrichment analysis e.g. 3702.

## target miRNA or gene
```
Rscript get_mirna_targets.R /path/to/miRNAlist.csv /path/to/output_folder
Rscript get_gene_targets.R --input gene_differents.csv --species mmu

```
