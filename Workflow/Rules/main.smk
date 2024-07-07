workdir: "/home/ubuntu/Desktop/Dataprocessing/"

configfile: "Config/config.yaml"
SAMPLES = config["samples"]
INDEX_PREFIX = "out/bowtie2_indices/genome"

include: "quality_check.smk"
include: "variant_calling.smk"
include: "downstream_processing.smk"


rule all:
    input:
        expand("out/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("out/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
        expand("out/samtools/mapped/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("out/varscan/{sample}.varScan.snp.filter", sample=SAMPLES),
        expand("out/varscan/{sample}.varScan.indel.filter", sample=SAMPLES),
        expand("out/varscan/{sample}.mpileup.readcounts", sample=SAMPLES),
        expand("out/mpileup_downstream/{sample}.var.flt.vcf", sample=SAMPLES),
        expand("out/samtools/{sample}.baq.var.raw.bcf", sample=SAMPLES),
        expand("out/images/QC_plot.png")
