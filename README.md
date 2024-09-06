## step1:
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

## step2: 

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







