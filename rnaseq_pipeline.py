#!/usr/bin/env python3

import os
import subprocess
import argparse
import glob

def run_star(fastq1, fastq2, index_dir, output_dir, threads, is_paired, sample_name):
    """ Run STAR to align FASTQ files to the reference genome. """

    # Create the STAR output directory if it doesn't exist
    sample_output_dir = os.path.join(output_dir, sample_name)
    if not os.path.exists(sample_output_dir):
        os.makedirs(sample_output_dir)

    # Prepare STAR command
    star_cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", index_dir,
        "--outFileNamePrefix", os.path.join(sample_output_dir, "STAR_"),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--quantMode", "GeneCounts"
    ]

    # Append FASTQ files based on whether it's paired-end or single-end
    if is_paired:
        star_cmd += ["--readFilesIn", fastq1, fastq2]
    else:
        star_cmd += ["--readFilesIn", fastq1]

    # Run STAR
    print(f"Running STAR for {sample_name} with command: {' '.join(star_cmd)}")
    subprocess.run(star_cmd, check=True)

    # Return the path to the generated BAM file
    return os.path.join(sample_output_dir, "STAR_Aligned.sortedByCoord.out.bam")


def run_featurecounts(bam_files, gtf_file, output_file, threads, is_paired):
    """ Run featureCounts to count reads mapped to features for multiple BAM files. """
    
    featurecounts_cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", gtf_file,
        "-o", output_file,
        "-t", "exon",
        "-g", "gene_id"
    ]

    if is_paired:
        featurecounts_cmd += ["-p"]

    # Add all BAM files
    featurecounts_cmd += bam_files

    # Run featureCounts
    print(f"Running featureCounts with command: {' '.join(featurecounts_cmd)}")
    subprocess.run(featurecounts_cmd, check=True)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="RNA-Seq Pipeline: Align multiple FASTQs and run featureCounts using STAR and featureCounts.")
    parser.add_argument("--fastq-dir", required=True, help="Directory containing FASTQ files.")
    parser.add_argument("--index-dir", required=True, help="Path to STAR genome index directory.")
    parser.add_argument("--gtf-file", required=True, help="Path to the GTF annotation file.")
    parser.add_argument("--output-dir", required=True, help="Directory to store the output files.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use.")
    
    args = parser.parse_args()

    # Check for paired-end or single-end FASTQ files
    fastq1_files = sorted(glob.glob(os.path.join(args.fastq_dir, "*_R1*.fastq*")))
    fastq2_files = sorted(glob.glob(os.path.join(args.fastq_dir, "*_R2*.fastq*")))

    # List to store BAM files
    bam_files = []

    # Loop through FASTQ files and run STAR for each pair
    for fastq1 in fastq1_files:
        sample_name = os.path.basename(fastq1).split("_R1")[0]  # Get the sample name
        fastq2 = None

        # Check if it's paired-end by looking for the corresponding R2 file
        if fastq2_files:
            fastq2 = fastq1.replace("_R1", "_R2")
            is_paired = os.path.exists(fastq2)
        else:
            is_paired = False

        # Run STAR alignment
        bam_file = run_star(fastq1, fastq2, args.index_dir, args.output_dir, args.threads, is_paired, sample_name)
        bam_files.append(bam_file)

    # The output file for featureCounts
    featurecounts_output = os.path.join(args.output_dir, "all_samples_gene_counts.txt")

    # Run featureCounts on all BAM files
    run_featurecounts(bam_files, args.gtf_file, featurecounts_output, args.threads, is_paired)


if __name__ == "__main__":
    main()
