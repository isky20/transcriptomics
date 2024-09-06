# rnaseq_pipeline.smk

# Load the config file
configfile: "config.yaml"

# Load samples from the CSV file
SAMPLES, FASTQ1, FASTQ2 = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

# Define the STAR rule
rule star_align:
    input:
        fastq1="data/fastq/{sample}_R1.fastq.gz",
        fastq2="data/fastq/{sample}_R2.fastq.gz"
    output:
        bam="results/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        counts="results/alignments/{sample}_ReadsPerGene.out.tab"
    params:
        index=config["index_dir"],
        prefix="results/alignments/{sample}_"
    threads: config["threads"]
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {params.index} \
             --readFilesIn {input.fastq1} {input.fastq2} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --outSAMunmapped Within \
             --outSAMattributes Standard
        """

# Define the featureCounts rule
rule featurecounts:
    input:
        bam="results/alignments/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        counts="results/counts/{sample}_gene_counts.txt"
    params:
        gtf=config["gtf_file"]
    threads: config["threads"]
    shell:
        """
        featureCounts -T {threads} \
                      -a {params.gtf} \
                      -o {output.counts} \
                      -t exon \
                      -g gene_id \
                      {input.bam}
        """

# Define the all rule to align all samples and run featureCounts
rule all:
    input:
        expand("results/alignments/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("results/counts/{sample}_gene_counts.txt", sample=SAMPLES)
