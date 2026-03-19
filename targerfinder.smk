#!/usr/bin/env python3

targetfinder = "docker://quay.io/biocontainers/targetfinder:1.7--0"

SCORE_CUTOFF  = 4.0
TARGET_DB     = "data/genome/JN3_transcript_clean.fa"  

grouped_srnas, = glob_wildcards("target_finder/grouped/{group}.fasta")
individual_srnas, = glob_wildcards("target_finder/individual/{sample}.fasta")

rule target:
    input:
        expand(
            "results/targetfinder/{sample}_targets.{fmt}",
            sample=grouped_srnas&individual_srnas,
            fmt=["table", "gff", "json"],
        ),
        "results/targetfinder/all_targets_combined.tsv",

rule targetfinder:
    input:
        srna_fa = "data/srnas/{sample}.fa",
        db      = TARGET_DB,
    output:
        table = "results/targetfinder/{sample}_targets.table",
        gff   = "results/targetfinder/{sample}_targets.gff",
        json  = "results/targetfinder/{sample}_targets.json",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 60,
    container:
        targetfinder
    shell:
        """
        SRNA_SEQ=$(grep -v '^>' {input.srna_fa} | tr -d '\n')
        SRNA_NAME="{wildcards.sample}"

        targetfinder.pl \
            -s "$SRNA_SEQ" \
            -d {input.db} \
            -q "$SRNA_NAME" \
            -c {SCORE_CUTOFF} \
            -p $FMT \
            -r \
            > {output}
        """

rule combine_results:
    input:
        expand("results/targetfinder/{sample}_targets.txt", sample=individual_srnas),
    output:
        "results/targetfinder/all_targets_combined.tsv",
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 15,
    run:
        written_header = False
        with open(output[0], "w") as out_fh:
            for f in sorted(input):
                with open(f) as in_fh:
                    lines = in_fh.readlines()
                if not lines:
                    continue
                if not written_header and lines[0].startswith("#"):
                    out_fh.write(lines[0])
                    written_header = True
                    lines = lines[1:]
                elif lines[0].startswith("#"):
                    lines = lines[1:]
                for line in lines:
                    out_fh.write(line)