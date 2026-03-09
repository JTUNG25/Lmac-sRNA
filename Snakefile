#!/usr/bin/env python3

fastp      = "docker://quay.io/biocontainers/fastp:1.0.1--heae3180_0"
bowtie1    = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"
shortstack = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"

input_samples, = glob_wildcards("raw_data_srna/links/{sample}_R1.fastq.gz")

rule target:
    input:
        expand("results/shortstack/{sample}/Results.txt", sample=input_samples)

rule run_shortstack:
    input:
        reads = "results/merged_srna/{sample}.merged.fastq.gz",
        index = rules.bowtie_index.output.idx
    output:
        "results/shortstack/{sample}/Results.txt"
    threads: 8
    container:
        shortstack
    shell:
        """
        ShortStack --readfile {input.reads} \
                   --bowtie_index data/bowtie_index \
                   --outdir results/shortstack/{wildcards.sample} \
                   --threads {threads} \
                   --mmap u
        """

rule bowtie_index:
    input:
        fasta = "data/JN3.fasta"
    output:
        idx = multiext("data/bowtie_index", ".1.ebwt", ".2.ebwt", ".3.ebwt",
                       ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
    container:
        bowtie1
    shell:
        "bowtie-build {input.fasta} data/bowtie_index"

rule fastp:
    input:
        r1 = "raw_data_srna/links/{sample}_R1.fastq.gz",
        r2 = "raw_data_srna/links/{sample}_R2.fastq.gz"
    output:
        merged = "results/merged_srna/{sample}.merged.fastq.gz",
        html   = "results/merged_srna/{sample}_fastp.html",
        json   = "results/merged_srna/{sample}_fastp.json"
    threads: 4
    container:
        fastp
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              --merge --merged_out {output.merged} \
              --detect_adapter_for_pe \
              --length_required 18 \
              --length_limit 30 \
              --thread {threads} \
              --html {output.html} --json {output.json}
        """