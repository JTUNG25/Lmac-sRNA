#!/usr/bin/env python3

fastp      = "docker://quay.io/biocontainers/fastp:1.0.1--heae3180_0"
bowtie1    = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"
shortstack = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"

wt_samples = [
    "D5-5-1", "D5-5-2", "D5-5-3",
    "D5-6-1", "D5-6-2", "D5-6-3",
    "D5-7-1", "D5-7-2", "D5-7-3",
]

rule target:
    input:
        "results/shortstack/wildtype_pooled/Results.txt"

rule bowtie_index:
    input:
        fasta = "data/genome/JN3.fasta"
    output:
        idx = multiext("data/bowtie_index", ".1.ebwt", ".2.ebwt", ".3.ebwt",
                       ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
    container:
        bowtie1
    shell:
        "bowtie-build {input.fasta} data/bowtie_index"

rule shortstack_pooled:
    input:
        reads = expand("data/fastp/merged/{sample}.merged.fastq.gz", sample=wt_samples),
        index = rules.bowtie_index.output.idx
    output:
        "results/shortstack/wildtype_pooled/Results.txt"
    threads: 8
    container:
        shortstack
    shell:
        """
        ShortStack --readfile {input.reads} \
                   --bowtie_index data/bowtie_index \
                   --outdir results/shortstack/wildtype_pooled \
                   --threads {threads} \
                   --mmap u
        """ 

rule fastp:
    input:
        r1 = "data/fastp/{sample}_R1.fastqsanger.gz",
        r2 = "data/fastp/{sample}_R2.fastqsanger.gz"
    output:
        merged = "data/fastp/merged/{sample}.merged.fastq.gz",
        html   = "data/fastp/merged/{sample}_fastp.html",
        json   = "data/fastp/merged/{sample}_fastp.json"
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