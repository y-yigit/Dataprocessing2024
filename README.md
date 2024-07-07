# Snakemake Fast-GBS

## About

Snakemake workflow for [Whole Exome Sequencing](https://bio-protocol.org/pdf/bio-protocol2805.pdf)

The snakemake pipeline was tested on SARS-CoV-2 strains which were downloaded from the SRA archive.
The fastq.gz files were extracted with: `decompress fastq seqtk seq {sample}.fastq.gz > {sample}.fastq`

The reference sequence for SARS-CoV-2 was downloaded from NCBI RefSeq.

## Requirements

- mamba 1.5.8
- conda 24.3.0
- Bash 5.0.17(1)-release
- Snakemake 8.14.0
- Python 3.12.3
- Bcftools 1.20
- FastQC v0.12.1
- Bowtie2 2.5.4
- Openjdk 22.0.
- Samtools 1.20
- VarScan.v2.3.9

## Notes

Setup environment: mamba create -c conda-forge -c bioconda -n snake_env snakemake snakemake-wrapper-utils fastqc bowtie2 samtools bcftools vcftools

Create DAG: `snakemake --snakefile Workflow/Rules/main.smk --dag | dot -Tpng > workflow_dag.png`

Currently, the rule visualise_varscan_output runs in incorrect order.
