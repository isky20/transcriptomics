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
Rscript deg_analysis.R --raw.reads.csv "featureCount.csv" \
                       --colname.case "CasC2_L3" \
                       --number.case.samples 4 \
                       --named.case "CasC" \
                       --colname.control "CasS2_L3" \
                       --number.control.samples 4 \
                       --named.control "CasS" \
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
Rscript deg_pd_analysis.R --count-data-file "data_pd.csv" \
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


## target miRNA
```
Rscript get_mirna_targets.R /path/to/miRNAlist.csv /path/to/output_folder
```
