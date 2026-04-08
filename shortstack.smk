#!/usr/bin/env python3

fastp = "docker://quay.io/biocontainers/fastp:1.0.1--heae3180_0"
shortstack = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"

wt_samples = [
    "D5-5-1",
    "D5-5-2",
    "D5-5-3",
    "D5-6-1",
    "D5-6-2",
    "D5-6-3",
    "D5-7-1",
    "D5-7-2",
    "D5-7-3",
]

(all_samples,) = glob_wildcards("raw_data_srna/links/{sample}_R1.fastq.gz")
all_samples = [s for s in all_samples if not s.startswith("RNA")]


rule target:
    input:
        expand("results/shortstack/mutants/{sample}/Results.txt", sample=all_samples),


rule shortstack_locifile:
    input:
        reads="data/merged_srna/{sample}.merged.fastq.gz",
        loci="results/shortstack/wt/Results.txt",
        genome="data/genome/JN3.fasta",
    output:
        "results/shortstack/mutants/{sample}/Results.txt",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    container:
        shortstack
    shell:
        """
        rm -rf results/shortstack/mutants/{wildcards.sample}
        ShortStack --readfile {input.reads} \
                   --genomefile {input.genome} \
                   --locifile {input.loci} \
                   --outdir results/shortstack/mutants/{wildcards.sample} \
                   --threads {threads} \
                   --mmap u
        """


rule shortstack_pooled:
    input:
        reads=expand("data/merged_srna/{sample}.merged.fastq.gz", sample=wt_samples),
        genome="data/genome/JN3.fasta",
    output:
        "results/shortstack/wt/Results.txt",
    threads: 16
    resources:
        mem_mb=32000,
        runtime=120,
    container:
        shortstack
    shell:
        """
        rm -rf results/shortstack/wt
        ShortStack --readfile {input.reads} \
                   --genomefile {input.genome} \
                   --outdir results/shortstack/wt \
                   --threads {threads} \
                   --dn_mirna \
                   --mmap u
        """


rule fastp:
    input:
        r1="raw_data_srna/links/{sample}_R1.fastq.gz",
        r2="raw_data_srna/links/{sample}_R2.fastq.gz",
    output:
        fastq_gz="data/merged_srna/{sample}.merged.fastq.gz",
        html="data/merged_srna/{sample}_fastp.html",
        json="data/merged_srna/{sample}_fastp.json",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    container:
        fastp
    log:
        "logs/fastp/{sample}.log",
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              --merge --merged_out {output.fastq_gz} \
              --detect_adapter_for_pe \
              --disable_quality_filtering \
              --length_required 18 \
              --length_limit 30 \
              --thread {threads} \
              --html {output.html} --json {output.json}
        """
