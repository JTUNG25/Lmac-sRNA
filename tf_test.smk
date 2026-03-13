#!/usr/bin/env python3

targetfinder = "docker://quay.io/biocontainers/targetfinder:1.7--0"

SCORE_CUTOFF = 3.0
TARGET_DB    = "data/genome/JN3_transcript_clean.fa"
TEST_SAMPLE  = "D1_top50"  

rule target:
    input:
        expand(
            "results/test/{sample}_targets.{fmt}",
            sample=[TEST_SAMPLE],
            fmt=["table", "gff", "json"],
        ),
        "results/test/hit_summary.txt",

rule targetfinder:
    input:
        srna_fa = "target_finder/individual/{sample}.fasta",
        db      = TARGET_DB,
    output:
        table = "results/test/{sample}_targets.table",
        gff   = "results/test/{sample}_targets.gff",
        json  = "results/test/{sample}_targets.json",
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

        for FMT in table gff json; do
            targetfinder.pl \
                -s "$SRNA_SEQ" \
                -d {input.db} \
                -q "$SRNA_NAME" \
                -c {SCORE_CUTOFF} \
                -p $FMT \
                > results/test/{wildcards.sample}_targets.$FMT
        done
        """

rule summarise_hits:
    input:
        table = "results/test/{sample}_targets.table".format(sample=TEST_SAMPLE),
    output:
        summary = "results/test/hit_summary.txt",
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 5,
    run:
        with open(input.table) as fh:
            lines = [l for l in fh if not l.startswith("#") and l.strip()]

        n_hits = len(lines)

        scores = []
        for line in lines:
            parts = line.split("\t")
            try:
                scores.append(float(parts[4]))
            except (IndexError, ValueError):
                pass

        with open(output.summary, "w") as out:
            out.write("=" * 50 + "\n")
            out.write("TargetFinder Test Run Summary\n")
            out.write("=" * 50 + "\n")
            out.write(f"Sample        : {TEST_SAMPLE}\n")
            out.write(f"Score cutoff  : {SCORE_CUTOFF}\n")
            out.write(f"Database      : {TARGET_DB}\n")
            out.write(f"Total hits    : {n_hits}\n")
            if scores:
                out.write(f"Score range   : {min(scores):.2f} - {max(scores):.2f}\n")
                out.write(f"Mean score    : {sum(scores)/len(scores):.2f}\n")
                out.write("\nScore distribution:\n")
                for threshold in [0.0, 1.0, 2.0, 3.0]:
                    count = sum(1 for s in scores if s <= threshold)
                    out.write(f"  <= {threshold:.1f}  :  {count} hits\n")
            out.write("\n--- Interpretation guide ---\n")
            out.write("0 hits          : cutoff too strict, try 4.0\n")
            out.write("1-20 hits       : good, high-confidence targets\n")
            out.write("20-100 hits     : reasonable, review score distribution\n")
            out.write(">100 hits       : consider tightening cutoff\n")
            out.write("=" * 50 + "\n")