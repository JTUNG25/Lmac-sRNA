#!/usr/bin/env python3

fastp = "docker://quay.io/biocontainers/fastp:1.0.1--heae3180_0"
shortstack = "docker://quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"
samtools = "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
bowtie = "docker://quay.io/biocontainers/bowtie:1.3.1--py310h4070885_6"

TRANSCRIPT_FASTA = "data/genome/JN3_transcript_clean.fa"

WT_SAMPLES = [
    "D5-5-1",
    "D5-5-2",
    "D5-5-3",
    "D5-6-1",
    "D5-6-2",
    "D5-6-3",
    "D5-7-1",
    "D5-7-2",
    "D5-7-3",
    "WT-1-1",
    "WT-1-2",
    "WT-1-3",
]


(all_samples,) = glob_wildcards("merged_raw/{sample}_R1.fastq.gz")
all_samples = [s for s in all_samples if not s.startswith("RNA")]


rule target:
    input:
        expand("results/shortstack/all/{sample}/Results.txt", sample=all_samples),
        expand("results/transcript_counts/{sample}.counts.txt", sample=all_samples),
        expand("results/transcript_counts/{sample}.21nt.counts.txt", sample=all_samples),
        expand("results/transcript_counts/{sample}.22nt.counts.txt", sample=all_samples),
        "results/transcript_counts/all_samples_counts.tsv",
        "results/transcript_counts/all_samples_21nt.tsv",
        "results/transcript_counts/all_samples_22nt.tsv",


rule shortstack_locifile:
    input:
        reads="data/fastp_srna/{sample}.merged.fastq.gz",
        loci="results/shortstack/wt/Results.txt",
        genome="data/genome/JN3.fasta",
    output:
        "results/shortstack/all/{sample}/Results.txt",
    container:
        shortstack
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        """
        rm -rf results/shortstack/all/{wildcards.sample}
        ShortStack --readfile {input.reads} \
                   --genomefile {input.genome} \
                   --locifile {input.loci} \
                   --outdir results/shortstack/all/{wildcards.sample} \
                   --threads {threads} \
                   --mmap u
        """


rule shortstack_pooled:
    input:
        reads=expand("data/fastp_srna/{sample}.merged.fastq.gz", sample=WT_SAMPLES),
        genome="data/genome/JN3.fasta",
    output:
        "results/shortstack/wt/Results.txt",
    container:
        shortstack
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
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


rule bowtie_index:
    input:
        TRANSCRIPT_FASTA,
    output:
        multiext(
            "data/genome/JN3_transcripts",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
    container:
        bowtie
    shell:
        """
        bowtie-build {input} data/genome/JN3_transcripts
        """


rule bowtie_map:
    input:
        reads="data/fastp_srna/{sample}.merged.fastq.gz",
        index=multiext(
            "data/genome/JN3_transcripts",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
    output:
        counts_total="results/transcript_counts/{sample}.counts.txt",
        counts_21="results/transcript_counts/{sample}.21nt.counts.txt",
        counts_22="results/transcript_counts/{sample}.22nt.counts.txt",
        stats="results/transcript_counts/{sample}.stats.txt",
    container:
        bowtie
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    shell:
        """
        bowtie -x data/genome/JN3_transcripts \
               -U {input.reads} \
               -p {threads} \
               --no-unal \
               -a --best --strata \
               --sam \
               2> {output.stats} \
        | awk '$2 != 4' \
        | tee \
            >(awk '{{if(length($10)==21) print $3}}' | sort | uniq -c | awk '{{print $2"\t"$1}}' > {output.counts_21}) \
            >(awk '{{if(length($10)==22) print $3}}' | sort | uniq -c | awk '{{print $2"\t"$1}}' > {output.counts_22}) \
        | awk '{{print $3}}' | sort | uniq -c | awk '{{print $2"\t"$1}}' \
        > {output.counts_total}
        """


rule merge_counts:
    input:
        total=expand(
            "results/transcript_counts/{sample}.counts.txt", sample=all_samples
        ),
        nt21=expand(
            "results/transcript_counts/{sample}.21nt.counts.txt", sample=all_samples
        ),
        nt22=expand(
            "results/transcript_counts/{sample}.22nt.counts.txt", sample=all_samples
        ),
    output:
        total="results/transcript_counts/all_samples_counts.tsv",
        nt21="results/transcript_counts/all_samples_21nt.tsv",
        nt22="results/transcript_counts/all_samples_22nt.tsv",
    run:
        import pandas as pd

        def merge_counts(files, samples):
            dfs = []
            for f, sample in zip(files, samples):
                df = pd.read_csv(
                    f, sep="\t", header=None, names=["transcript", sample]
                )
                dfs.append(df.set_index("transcript"))
            merged = pd.concat(dfs, axis=1).fillna(0).astype(int)
            merged.index.name = "transcript"
            return merged

        merge_counts(input.total, all_samples).to_csv(output.total, sep="\t")
        merge_counts(input.nt21, all_samples).to_csv(output.nt21, sep="\t")
        merge_counts(input.nt22, all_samples).to_csv(output.nt22, sep="\t")


rule fastp:
    input:
        r1="merged_raw/{sample}_R1.fastq.gz",
        r2="merged_raw/{sample}_R2.fastq.gz",
    output:
        fastq_gz="data/fastp_srna/{sample}.merged.fastq.gz",
        html="data/fastp_srna/{sample}_fastp.html",
        json="data/fastp_srna/{sample}_fastp.json",
    log:
        "logs/fastp/{sample}.log",
    container:
        fastp
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
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
