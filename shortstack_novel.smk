#!/usr/bin/env python3

fastp      = "docker://quay.io/biocontainers/fastp:1.0.1--heae3180_0"
shortstack = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"


(all_samples,) = glob_wildcards("raw_data_srna_2/merged/{sample}_R1.fastq.gz")

# WT controls are included so novel clusters are always called in context of the WT baseline 

GROUPS = {
    "dcl": {
        "samples": ["D1", "D2-2", "D2-3", "D2-4",
                    "D5-5", "D5-6"],          # WT controls for batch1 + batch2
    },
    "ago": {
        "samples": ["A1-1", "A1-2", "A1-3", "A3",
                    "A13-1", "A13-2",
                    "D5-6", "WT-1"],          # WT controls for batch2 + batch4
    },
    "rdrp": {
        "samples": ["R1", "R2-2", "R2-3", "R2-4", "R2-5",
                    "R12-1", "R12-2", "R12-3",
                    "D5-7", "WT-1"],          # WT controls for batch3 + batch4
    },
}


rule target:
    input:
        expand("results/shortstack/novel/{group}/Results.txt",
               group=GROUPS.keys()),


rule shortstack_novel:
    input:
        reads  = lambda wc: expand(
                     "data/fastp_srna/{sample}.merged.fastq.gz",
                     sample=GROUPS[wc.group]["samples"]
                 ),
        genome = "data/genome/JN3.fasta",
    output:
        "results/shortstack/novel/{group}/Results.txt",
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
    container:
        shortstack
    shell:
        """
        rm -rf results/shortstack/novel/{wildcards.group}
        ShortStack --readfile {input.reads} \
                   --genomefile {input.genome} \
                   --outdir results/shortstack/novel/{wildcards.group} \
                   --threads {threads} \
                   --dn_mirna \
                   --mmap u
        """


rule fastp:
    input:
        r1 = "merged_raw/{sample}_R1.fastq.gz",
        r2 = "merged_raw/{sample}_R2.fastq.gz",
    output:
        fastq_gz = "data/fastp_srna/{sample}.merged.fastq.gz",
        html     = "data/fastp_srna/{sample}_fastp.html",
        json     = "data/fastp_srna/{sample}_fastp.json",
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
