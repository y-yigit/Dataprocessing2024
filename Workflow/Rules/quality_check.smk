# Preprocessing of raw data by samtools and bowtie
# Preprocessing of raw data by samtools and bowtie
configfile: "Config/config.yaml"
OUTPUT_DIR = config["output_dir"]
REF_GENOME = config["reference"]
INDEX_PREFIX = "out/bowtie2_indices/genome"

# Rule to run FastQC for quality control of raw sequencing reads
rule fastqc:
    input:
        "data/{sample}.fastq",
    output:
        html="out/fastqc/{sample}_fastqc.html",
        zip = "out/fastqc/{sample}_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc -o {OUTPUT_DIR}/fastqc {input}"

# Rule to build a Bowtie2 index from a reference genome
rule bowtie2_build:
    input:
        reference = "data/reference.fasta",
    output:
        INDEX_PREFIX + ".1.bt2",
        INDEX_PREFIX + ".2.bt2",
        INDEX_PREFIX + ".3.bt2",
        INDEX_PREFIX + ".4.bt2",
        INDEX_PREFIX + ".rev.1.bt2",
        INDEX_PREFIX + ".rev.2.bt2",
        #"out/bowtie2_alignments",
    log:
        "logs/bowtie2/build.log"
    shell:
        """
        bowtie2-build {input.reference} {INDEX_PREFIX}
        mkdir -p out/bowtie2_alignments
        """

# Rule to align sequencing reads to a reference genome using Bowtie2
rule bowtie2:
    input:
        reads = "data/{sample}.fastq",
        index_files = [
            INDEX_PREFIX + ".1.bt2",
            INDEX_PREFIX + ".2.bt2",
            INDEX_PREFIX + ".3.bt2",
            INDEX_PREFIX + ".4.bt2",
            INDEX_PREFIX + ".rev.1.bt2",
            INDEX_PREFIX + ".rev.2.bt2"
        ]
    output:
        "out/bowtie2_alignments/{sample}_aligned_reads.sam"
    log:
        "logs/bowtie2/alignment_{sample}.log"
    shell:
        """
        bowtie2 -x {INDEX_PREFIX} -U {input.reads} -S {output}
        """

rule sam_to_bam:
    input:
        "out/bowtie2_alignments/{sample}_aligned_reads.sam"
    output:
        "out/bowtie2_alignments/{sample}_aligned_reads.bam"
    log:
        "logs/bowtie2/{sample}_to_bam.log"
    shell:
        "samtools view -bS {input} > {output}"

# Rule to sort BAM files using samtools
rule samtools_sort:
    input:
        "out/bowtie2_alignments/{sample}_aligned_reads.bam"
    output:
        "out/samtools/mapped/{sample}.sorted.bam",
    log:
        "logs/samtools/mapped_{sample}.log"
    shell:
        """
        samtools sort -o {output} {input}
        """

# Rule to index sorted BAM files using samtools
rule samtools_index:
    input:
        "out/samtools/mapped/{sample}.sorted.bam",
    output:
        "out/samtools/mapped/{sample}.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    shell:
        "samtools index {input} -o {output}"

# Rule to generate a pileup file using samtools mpileup
rule mpilup:
    input:
        # single or list of bam files
        "out/samtools/mapped/{sample}.sorted.bam",
    output:
        "out/mpileup/{sample}.mpileup",
    log:
        "logs/samtools_mpileup/{sample}.log",
    shell:
        "mkdir -p out/mpileup && samtools mpileup -f data/{REF_GENOME}.fasta {input} > {output}"
