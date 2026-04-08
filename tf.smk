#!/usr/bin/env python3

targetfinder = "docker://quay.io/biocontainers/targetfinder:1.7--3"

SCORE_CUTOFF = 4.0
TARGET_DB = "data/genome/JN3_transcript_clean.fa"

TEST_SAMPLES = [
    "D1_top50",
    "D2-2_top50",
    "D2-3_top50",
    "D2-4_top50",
    "A1-1_top50",
    "A1-2_top50",
    "A1-3_top50",
    "A3_top50",
    "R1_top50",
    "R2-2_top50",
    "R2-3_top50",
    "R2-4_top50",
    "R2-5_top50",
]


rule target:
    input:
        expand(
            "results/tf_edger/{sample}_targets.{fmt}",
            sample=TEST_SAMPLES,
            fmt=["table", "gff", "json"],
        ),
        "results/tf_edger/hit_summary.txt",


rule targetfinder:
    input:
        srna_fa="data/edger_seq/individual/{sample}.fasta",
        db=TARGET_DB,
    output:
        table="results/tf_edger/{sample}_targets.table",
        gff="results/tf_edger/{sample}_targets.gff",
        json="results/tf_edger/{sample}_targets.json",
    container:
        targetfinder
    threads: 1
    resources:
        mem_mb=32000,
        runtime=240,
    params:
        score=SCORE_CUTOFF,
    shell:
        """
        SRNA_NAME=""
        SRNA_SEQ=""

        while IFS= read -r line || [[ -n "$line" ]]; do
            if [[ "$line" == ">"* ]]; then
                if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
                    for FMT in table gff json; do
                        targetfinder.pl \
                            -s "$SRNA_SEQ" \
                            -d {input.db} \
                            -q "$SRNA_NAME" \
                            -c {params.score} \
                            -p $FMT \
                            -t 1 -r \
                            >> results/tf_edger/{wildcards.sample}_targets.$FMT
                    done
                fi
                SRNA_NAME="${{line:1}}"
                SRNA_SEQ=""
            else
                SRNA_SEQ="${{SRNA_SEQ}}${{line}}"
            fi
        done < {input.srna_fa}

        # Flush the final sequence
        if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
            for FMT in table gff json; do
                targetfinder.pl \
                    -s "$SRNA_SEQ" \
                    -d {input.db} \
                    -q "$SRNA_NAME" \
                    -c {params.score} \
                    -p $FMT \
                    -t 1 -r \
                    >> results/tf_edger/{wildcards.sample}_targets.$FMT
            done
        fi
        """


rule summarise_hits:
    input:
        tables=expand("results/tf_edger/{sample}_targets.table", sample=TEST_SAMPLES),
    output:
        summary="results/tf_edger/hit_summary.txt",
    run:
        sample_results = {}

        for table_file in input.tables:
            sample_name = table_file.split("/")[-1].replace("_targets.table", "")
            with open(table_file) as f:
                data_lines = [l for l in f if not l.startswith("#") and l.strip()]
            scores = []
            for line in data_lines:
                parts = line.split("\t")
                try:
                    scores.append(float(parts[4]))
                except (IndexError, ValueError):
                    pass
            sample_results[sample_name] = scores

        total_hits = sum(len(s) for s in sample_results.values())

        with open(output.summary, "w") as out:
            out.write("TargetFinder Results Summary\n")
            out.write("=" * 40 + "\n")
            out.write(f"Database    : {TARGET_DB}\n")
            out.write(f"Score cutoff: {SCORE_CUTOFF}\n")
            out.write(f"Samples     : {len(TEST_SAMPLES)}\n\n")

            for sample_name, scores in sample_results.items():
                hits = len(scores)
                if scores:
                    out.write(
                        f"{sample_name}: {hits} hits  (scores {min(scores):.1f} – {max(scores):.1f})\n"
                    )
                else:
                    out.write(f"{sample_name}: 0 hits\n")

            out.write(f"\nTotal hits          : {total_hits}\n")
            out.write(
                f"Average per sample  : {total_hits/ len(TEST_SAMPLES):.1f}\n"
            )

            no_hits = sum(1 for s in sample_results.values() if len(s) == 0)
            few_hits = sum(1 for s in sample_results.values() if 1 <= len(s) <= 10)
            many_hits = sum(1 for s in sample_results.values() if len(s) > 10)

            out.write(f"\nSample distribution:\n")
            out.write(f"  No hits  : {no_hits}\n")
            out.write(f"  1–10 hits: {few_hits}\n")
            out.write(f"  >10 hits : {many_hits}\n")
