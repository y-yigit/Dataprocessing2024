# Variant calling by mpileup and varscan

# Rule to call SNPs using VarScan
rule varscan_snp:
    input:
        "out/mpileup/{sample}.mpileup",
    output:
        "out/varscan/{sample}.varScan.snp",
    log:
        "logs/varscan_mpileup_snp/{sample}.log",
    shell:
        "mkdir -p out/varscan && java -jar tools/VarScan.v2.3.9.jar mpileup2snp {input} > {output}"

# Rule to call indels using VarScan
rule varscan_indel:
    input:
        "out/mpileup/{sample}.mpileup",
    output:
        "out/varscan/{sample}.varScan.indel",
    log:
        "logs/varscan_mpileup_indel/{sample}.log",
    shell:
        "java -jar tools/VarScan.v2.3.9.jar mpileup2indel {input} > {output}"

# The rules below can produce empty output

# Rule to filter SNP calls using VarScan
rule varscan_snp_filter:
    input:
        snp="out/varscan/{sample}.varScan.snp",
        indel="out/varscan/{sample}.varScan.indel",
    output:
        "out/varscan/{sample}.varScan.snp.filter",
    log:
        "logs/varscan_snp_filter/{sample}.log",
    shell:
        "java -jar tools/VarScan.v2.3.9.jar filter {input.snp} --indel-file {input.indel} --output-file {output}"

# Rule to filter indel calls using VarScan
rule varscan_indel_filter:
    input:
       "out/varscan/{sample}.varScan.indel",
    output:
        "out/varscan/{sample}.varScan.indel.filter",
    log:
        "logs/varscan_indel_filter/{sample}.log",
    shell:
        "java -jar tools/VarScan.v2.3.9.jar filter {input} --output-file {output}"

# Rule to generate read counts using VarScan
rule varscan_read_counts:
    input:
        "out/mpileup/{sample}.mpileup", #.sam?
    output:
        "out/varscan/{sample}.mpileup.readcounts",
    log:
        "logs/varscan_readcounts/{sample}.log",
    shell:
        "java -jar tools/VarScan.v2.3.9.jar readcounts {input} > {output}"

rule visualise_varscan_output:
    output:
        "out/images/QC_plot.png"
    shell:
        "mkdir -p out/images && python3 Workflow/Python/QC.py"

