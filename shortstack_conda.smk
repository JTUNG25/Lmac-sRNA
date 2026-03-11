#!/usr/bin/env python3

wt_samples = [
    "D5-5-1", "D5-5-2", "D5-5-3",
    "D5-6-1", "D5-6-2", "D5-6-3",
    "D5-7-1", "D5-7-2", "D5-7-3",
]

rule target:
    input:
        "results/shortstack/de_novo/Results.txt"

rule fastp:
    input:
        r1 = "data/fastp/{sample}_R1.fastqsanger.gz",
        r2 = "data/fastp/{sample}_R2.fastqsanger.gz"
    output:
        merged = "data/merged/{sample}.merged.fastq.gz",
        html   = "data/merged/{sample}_fastp.html",
        json   = "data/merged/{sample}_fastp.json"
    threads: 4
    conda: "envs/fastp.yaml"
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

rule shortstack_pooled:
    input:
        reads = expand("data/merged/{sample}.merged.fastq.gz", sample=wt_samples),
        genome = "data/genome/JN3.fasta"
    output:
        "results/shortstack/de_novo/Results.txt"
    threads: 8
    conda: "envs/st.yaml"
    shell:
        """
        rm -rf results/shortstack/de_novo/
        ShortStack --readfile {input.reads} \
                   --genomefile {input.genome} \
                   --outdir results/shortstack/de_novo \
                   --threads {threads} \
                   --dn_mirna \
                   --mmap u
        """