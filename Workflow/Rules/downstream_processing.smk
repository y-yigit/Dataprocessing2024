# Downstream processing of files by samtools and bcftoolss
configfile: "Config/config.yaml"

# Rule to generate a pileup file using samtools mpileup and bcftools

rule mpilup_bcf:
    input:
        reference="data/reference.fasta",
        bam="out/samtools/mapped/{sample}.sorted.bam",
    output:
        "out/mpileup_downstream/{sample}.var.raw.bcf",
    log:
        "logs/bcf_mpileup/{sample}.log",
    shell:
        "bcftools mpileup -f {input.reference} {input.bam} > {output}"

# Rule to filter variants using bcftools and vcfutils.pl
rule bcftools:
    input:
        "out/mpileup_downstream/{sample}.var.raw.bcf",
    output:
        "out/mpileup_downstream/{sample}.var.flt.vcf",
    log:
        "logs/vcf_mpileup/{sample}.log",
    shell:
        "bcftools view {input} | vcfutils.pl varFilter -D100 > {output}"

rule samtools_calmd:
    input:
        bam="out/samtools/mapped/{sample}.sorted.bam",
        reference="data/reference.fasta",
    output:
        "out/samtools/{sample}.baq.bam",
    log:
        "logs/samtools_calmd/{sample}.log",
    shell:
        "samtools calmd -Abr {input.bam} {input.reference} > {output}"

# Rule to generate a pileup file using samtools mpileup and bcftools after BAQ calculation
rule mpilup_downstream_bam:
    input:
        reference="data/reference.fasta",
        bam="out/samtools/{sample}.baq.bam",
    output:
        "out/samtools/{sample}.baq.var.raw.bcf",
    log:
        "logs/samtools/mpileup/{sample}.log",
    shell:
        "bcftools mpileup -f {input.reference} {input.bam} > {output}"
