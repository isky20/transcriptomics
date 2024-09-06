## step1:
```
python rnaseq_pipeline.py --fastq1 path/to/read1.fastq \
                          --fastq2 path/to/read2.fastq (if paired-end) \
                          --index-dir path/to/star/index \
                          --gtf-file path/to/annotation.gtf \
                          --output-dir path/to/output_directory \
                          --threads 8
```
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








