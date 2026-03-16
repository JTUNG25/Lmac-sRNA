#!/usr/bin/env python3

targetfinder = "docker://quay.io/biocontainers/targetfinder:1.7--0"

SCORE_CUTOFF = 5.0
TARGET_DB    = "data/genome/JN3_transcript_clean.fa"
TEST_SAMPLES = ["D1_top50", "D2-2_top50", "D2-3_top50", "D2-4_top50",
"A1-1_top50", "A1-2_top50", "A1-3_top50", "A3_top50",
"R1_top50", "R2-2_top50", "R2-3_top50", "R2-4_top50", "R2-5_top50"]


rule target:
    input:
        expand(
            "results/test/{sample}_targets.{fmt}",
            sample=TEST_SAMPLES,
            fmt=["table", "gff", "json"],
        ),
        expand("results/test/{sample}_diagnostics.txt", sample=TEST_SAMPLES),
        "results/test/hit_summary.txt",
        "results/test/workflow_diagnostics.txt",


rule targetfinder:
    input:
        srna_fa = "target_finder/individual/{sample}.fasta",
        db      = TARGET_DB,
    output:
        table = "results/test/{sample}_targets.table",
        gff   = "results/test/{sample}_targets.gff",
        json  = "results/test/{sample}_targets.json",
    params:
        score = SCORE_CUTOFF,
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 60,
    container:
        targetfinder
    log:
        "logs/targetfinder/{sample}.log",
    shell:
        """
        echo "========================================" > {log}
        echo "TargetFinder Analysis: {wildcards.sample}" >> {log}
        echo "========================================" >> {log}
        echo "Input FASTA: {input.srna_fa}" >> {log}
        echo "Database: {input.db}" >> {log}
        echo "Score cutoff: {params.score}" >> {log}
        echo "Start time: $(date)" >> {log}
        echo "" >> {log}

        # Check input file
        echo "Input file check:" >> {log}
        wc -l {input.srna_fa} >> {log}
        head -4 {input.srna_fa} >> {log}
        echo "" >> {log}

        SRNA_COUNT=$(grep -c "^>" {input.srna_fa})
        echo "Number of sRNA sequences: $SRNA_COUNT" >> {log}
        echo "" >> {log}

        # Initialise empty output files
        for FMT in table gff json; do
            > results/test/{wildcards.sample}_targets.$FMT
        done

        # Process one sequence at a time from the multi-FASTA
        SRNA_NAME=""
        SRNA_SEQ=""

        while IFS= read -r line || [[ -n "$line" ]]; do
            if [[ "$line" == ">"* ]]; then
                # Flush the previous sequence through TargetFinder
                if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
                    echo "Processing: $SRNA_NAME" >> {log}
                    for FMT in table gff json; do
                        targetfinder.pl \
                            -s "$SRNA_SEQ" \
                            -d {input.db} \
                            -q "$SRNA_NAME" \
                            -c {params.score} \
                            -p $FMT \
                            >> results/test/{wildcards.sample}_targets.$FMT \
                            2>> {log}
                    done
                fi
                # Start buffering the next sequence
                SRNA_NAME="${{line:1}}"
                SRNA_SEQ=""
            else
                SRNA_SEQ="${{SRNA_SEQ}}${{line}}"
            fi
        done < {input.srna_fa}

        # Flush the final sequence (not caught by the loop above)
        if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
            echo "Processing: $SRNA_NAME" >> {log}
            for FMT in table gff json; do
                targetfinder.pl \
                    -s "$SRNA_SEQ" \
                    -d {input.db} \
                    -q "$SRNA_NAME" \
                    -c {params.score} \
                    -p $FMT \
                    >> results/test/{wildcards.sample}_targets.$FMT \
                    2>> {log}
            done
        fi

        echo "" >> {log}
        echo "End time: $(date)" >> {log}
        echo "Table hits: $(grep -vc '^#' results/test/{wildcards.sample}_targets.table || echo 0)" >> {log}
        echo "Status: SUCCESS" >> {log}
        """


rule diagnose_sample:
    input:
        srna_fa = "target_finder/individual/{sample}.fasta",
        table   = "results/test/{sample}_targets.table",
        log     = "logs/targetfinder/{sample}.log",
    output:
        diagnostics = "results/test/{sample}_diagnostics.txt",
    run:
        with open(output.diagnostics, "w") as out:
            out.write("=" * 60 + "\n")
            out.write(f"DIAGNOSTIC REPORT: {wildcards.sample}\n")
            out.write("=" * 60 + "\n\n")

            with open(input.srna_fa) as f:
                srna_lines = f.readlines()
                srna_count = sum(1 for l in srna_lines if l.startswith(">"))
                seq_count  = sum(1 for l in srna_lines if not l.startswith(">") and l.strip())

            out.write("INPUT STATISTICS:\n")
            out.write("-" * 60 + "\n")
            out.write(f"Input file      : {input.srna_fa}\n")
            out.write(f"sRNA sequences  : {srna_count}\n")
            out.write(f"Sequence lines  : {seq_count}\n")
            out.write(f"Total file lines: {len(srna_lines)}\n\n")

            with open(input.table) as f:
                table_lines = f.readlines()

            data_lines = [l for l in table_lines if not l.startswith("#") and l.strip()]

            out.write("OUTPUT STATISTICS:\n")
            out.write("-" * 60 + "\n")
            out.write(f"Total table lines      : {len(table_lines)}\n")
            out.write(f"Data lines (hits)      : {len(data_lines)}\n\n")

            scores  = []
            targets = []
            for line in data_lines:
                parts = line.split("\t")
                try:
                    scores.append(float(parts[4]))
                    targets.append(parts)
                except (IndexError, ValueError):
                    pass

            out.write("PREDICTION STATISTICS:\n")
            out.write("-" * 60 + "\n")

            if scores:
                out.write(f"Total predictions: {len(scores)}\n")
                out.write(f"Score range      : {min(scores):.2f} - {max(scores):.2f}\n")
                out.write(f"Mean score       : {sum(scores)/len(scores):.2f}\n")
                out.write(f"Median score     : {sorted(scores)[len(scores)//2]:.2f}\n\n")
                out.write("SCORE DISTRIBUTION:\n")
                for threshold in [1.0, 2.0, 3.0, 4.0]:
                    count = sum(1 for s in scores if s <= threshold)
                    pct   = (count / len(scores)) * 100
                    out.write(f"  <= {threshold:.1f}: {count:4d} predictions ({pct:5.1f}%)\n")
                out.write("\nTOP 5 PREDICTIONS:\n")
                for i, (score, target) in enumerate(
                    sorted(zip(scores, targets), key=lambda x: x[0])[:5]
                ):
                    out.write(f"  {i+1}. Score: {score:.2f} - Target: {target[0] if target else 'N/A'}\n")
            else:
                out.write("NO PREDICTIONS FOUND\n")
                out.write("  - Score cutoff may be too strict\n")
                out.write("  - sRNAs may not have good matches in database\n")
                out.write("  - Try increasing cutoff to 5.0\n")

            out.write("\n" + "=" * 60 + "\n")
            out.write("INTERPRETATION GUIDE:\n")
            out.write("=" * 60 + "\n")
            out.write("0 hits      : Cutoff too strict, try 5.0\n")
            out.write("1-10 hits   : High-confidence targets\n")
            out.write("10-50 hits  : Good, selective targeting\n")
            out.write("50-100 hits : Reasonable, check score distribution\n")
            out.write(">100 hits   : Consider tightening cutoff\n")


rule summarise_hits:
    input:
        tables = expand("results/test/{sample}_targets.table", sample=TEST_SAMPLES),
        diags  = expand("results/test/{sample}_diagnostics.txt", sample=TEST_SAMPLES),
    output:
        summary = "results/test/hit_summary.txt",
    run:
        with open(output.summary, "w") as out:
            out.write("=" * 60 + "\n")
            out.write("TargetFinder Test Run - Overall Summary\n")
            out.write("=" * 60 + "\n\n")
            out.write(f"Score cutoff  : {SCORE_CUTOFF}\n")
            out.write(f"Database      : {TARGET_DB}\n")
            out.write(f"Samples tested: {', '.join(TEST_SAMPLES)}\n\n")

            total_hits = 0
            all_scores = []

            for table_file in input.tables:
                sample_name = table_file.split("/")[-1].replace("_targets.table", "")
                with open(table_file) as f:
                    lines = [l for l in f if not l.startswith("#") and l.strip()]
                scores = []
                for line in lines:
                    parts = line.split("\t")
                    try:
                        scores.append(float(parts[4]))
                        all_scores.append(float(parts[4]))
                    except (IndexError, ValueError):
                        pass
                total_hits += len(scores)
                out.write(f"{sample_name}:\n")
                out.write(f"  Hits: {len(scores)}\n")
                if scores:
                    out.write(f"  Score range : {min(scores):.2f} - {max(scores):.2f}\n")
                    out.write(f"  Mean score  : {sum(scores)/len(scores):.2f}\n")
                    for threshold in [1.0, 2.0, 3.0, 4.0]:
                        count = sum(1 for s in scores if s <= threshold)
                        out.write(f"  <= {threshold:.1f}      : {count} hits\n")
                else:
                    out.write("  No targets found\n")
                out.write("\n")

            out.write(f"Total hits across all samples: {total_hits}\n")
            if all_scores:
                out.write(f"Overall score range         : {min(all_scores):.2f} - {max(all_scores):.2f}\n")
                out.write(f"Overall mean score          : {sum(all_scores)/len(all_scores):.2f}\n")
            out.write("\n")
            out.write("Interpretation guide:\n")
            out.write("  0 hits      : cutoff too strict, try 5.0\n")
            out.write("  1-10 hits   : high confidence targets\n")
            out.write("  10-50 hits  : good, selective targeting\n")
            out.write("  50-200 hits : reasonable, check score distribution\n")
            out.write("  >200 hits   : consider tightening cutoff\n")


rule workflow_diagnostics:
    input:
        tables = expand("results/test/{sample}_targets.table", sample=TEST_SAMPLES),
        logs   = expand("logs/targetfinder/{sample}.log", sample=TEST_SAMPLES),
    output:
        wf_diag = "results/test/workflow_diagnostics.txt",
    run:
        with open(output.wf_diag, "w") as out:
            out.write("=" * 70 + "\n")
            out.write("WORKFLOW DIAGNOSTIC REPORT\n")
            out.write("=" * 70 + "\n\n")
            out.write(f"Timestamp   : {__import__('datetime').datetime.now()}\n")
            out.write(f"Container   : {targetfinder}\n")
            out.write(f"Score cutoff: {SCORE_CUTOFF}\n")
            out.write(f"Database    : {TARGET_DB}\n")
            out.write(f"Samples     : {len(TEST_SAMPLES)}\n\n")

            out.write("LOG FILE SUMMARIES:\n")
            out.write("-" * 70 + "\n")
            for log_file, sample in zip(input.logs, TEST_SAMPLES):
                out.write(f"\n{sample}:\n")
                try:
                    with open(log_file) as f:
                        log_content = f.read()
                    if "SUCCESS" in log_content:
                        out.write("  Status: SUCCESS\n")
                    elif "ERROR" in log_content:
                        out.write("  Status: ERROR\n")
                    else:
                        out.write("  Status: UNKNOWN\n")
                    for line in log_content.split("\n"):
                        if "sRNA sequences:" in line or "Table hits:" in line:
                            out.write(f"  {line.strip()}\n")
                except FileNotFoundError:
                    out.write("  Log file not found\n")

            out.write("\n" + "=" * 70 + "\n")
            out.write("RECOMMENDATIONS:\n")
            out.write("=" * 70 + "\n")
            out.write("1. Check individual sample diagnostic files\n")
            out.write("2. Review hit_summary.txt for overall statistics\n")
            out.write("3. Adjust score cutoff if needed\n")
            out.write("4. If satisfied, run full workflow on all samples\n")
            out.write("=" * 70 + "\n")