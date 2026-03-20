#!/usr/bin/env python3

targetfinder = "docker://quay.io/biocontainers/targetfinder:1.7--0"

SCORE_CUTOFF = 4.0
TARGET_DB = "data/genome/JN3.fasta"
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
        expand("results/tf/{sample}_targets.{fmt}", sample=TEST_SAMPLES, fmt=["table", "gff", "json"]),
        expand("results/tf/{sample}_diagnostics.txt", sample=TEST_SAMPLES),
        "results/tf/hit_summary.txt",
        "results/tf/workflow_diagnostics.txt",

rule targetfinder:
    input:
        srna_fa="data/shortstack_seq/individual/{sample}.fasta",
        db=TARGET_DB,
    output:
        table="results/tf/{sample}_targets.table",
        gff="results/tf/{sample}_targets.gff", 
        json="results/tf/{sample}_targets.json",
    params:
        score=SCORE_CUTOFF,
    threads: 1
    resources:
        mem_mb=32000,
        runtime=240,
    container: targetfinder
    log: "logs/targetfinder/{sample}.log",
    shell:
        """
        echo "TargetFinder Analysis: {wildcards.sample}" > {log}
        echo "=======================================" >> {log}

        # Find the correct TargetFinder executable
        TARGETFINDER_CMD=""
        for cmd in "targetfinder.pl" "targetfinder" "TargetFinder" "targetfinder-linux"; do
            if command -v "$cmd" > /dev/null 2>&1; then
                TARGETFINDER_CMD="$cmd"
                echo "Found TargetFinder executable: $cmd" >> {log}
                break
            fi
        done
        
        if [[ -z "$TARGETFINDER_CMD" ]]; then
            echo "ERROR: No TargetFinder executable found!" >> {log}
            ls /usr/local/bin/ >> {log} 2>&1 || echo "Cannot list /usr/local/bin/" >> {log}
            exit 1
        fi

        # Initialize output files
        for FMT in table gff json; do
            > results/tf/{wildcards.sample}_targets.$FMT
        done

        # Process sequences
        while IFS= read -r line || [[ -n "$line" ]]; do
            if [[ "$line" == ">"* ]]; then
                if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
                    echo "Processing: $SRNA_NAME" >> {log}
                    for FMT in table gff json; do
                        $TARGETFINDER_CMD -s "$SRNA_SEQ" -d {input.db} -q "$SRNA_NAME" -c {params.score} -p $FMT -t 1 -f -u 25 -r >> results/tf/{wildcards.sample}_targets.$FMT 2>> {log}
                    done
                fi
                SRNA_NAME="${{line:1}}"
                SRNA_SEQ=""
            else
                SRNA_SEQ="${{SRNA_SEQ}}${{line}}"
            fi
        done < {input.srna_fa}

        # Process final sequence
        if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
            echo "Processing: $SRNA_NAME" >> {log}
            for FMT in table gff json; do
                $TARGETFINDER_CMD -s "$SRNA_SEQ" -d {input.db} -q "$SRNA_NAME" -c {params.score} -p $FMT -t 1 -f -u 25 -r >> results/tf/{wildcards.sample}_targets.$FMT 2>> {log}
            done
        fi

        echo "Analysis complete" >> {log}
        """

rule diagnose_sample:
    input:
        srna_fa="data/shortstack_seq/individual/{sample}.fasta",
        table="results/tf/{sample}_targets.table",
        log="logs/targetfinder/{sample}.log",
    output:
        diagnostics="results/tf/{sample}_diagnostics.txt",
    run:
        with open(output.diagnostics, "w") as out:
            out.write(f"DIAGNOSTIC REPORT: {wildcards.sample}\n")
            out.write("=" * 40 + "\n")

            with open(input.srna_fa) as f:
                srna_count = sum(1 for l in f if l.startswith(">"))

            with open(input.table) as f:
                data_lines = [l for l in f if not l.startswith("#") and l.strip()]

            out.write(f"Input sequences: {srna_count}\n")
            out.write(f"Output hits: {len(data_lines)}\n")

            if data_lines:
                scores = []
                for line in data_lines:
                    parts = line.split("\t")
                    try:
                        scores.append(float(parts[4]))
                    except (IndexError, ValueError):
                        pass
                
                if scores:
                    out.write(f"Score range: {min(scores):.2f} - {max(scores):.2f}\n")
            else:
                out.write("NO HITS FOUND - try score cutoff 5.0 or 6.0\n")

rule summarise_hits:
    input:
        tables=expand("results/tf/{sample}_targets.table", sample=TEST_SAMPLES),
    output:
        summary="results/tf/hit_summary.txt",
    run:
        total_hits = 0
        with open(output.summary, "w") as out:
            out.write("TargetFinder Summary\n")
            out.write("=" * 20 + "\n")
            
            for table_file in input.tables:
                sample_name = table_file.split("/")[-1].replace("_targets.table", "")
                with open(table_file) as f:
                    hits = len([l for l in f if not l.startswith("#") and l.strip()])
                total_hits += hits
                out.write(f"{sample_name}: {hits} hits\n")
            
            out.write(f"Total hits: {total_hits}\n")

rule workflow_diagnostics:
    input:
        logs=expand("logs/targetfinder/{sample}.log", sample=TEST_SAMPLES),
    output:
        wf_diag="results/tf/workflow_diagnostics.txt",
    run:
        with open(output.wf_diag, "w") as out:
            out.write("Workflow Status\n")
            out.write("=" * 15 + "\n")
            
            for log_file, sample in zip(input.logs, TEST_SAMPLES):
                try:
                    with open(log_file) as f:
                        content = f.read()
                        if "Found TargetFinder executable" in content:
                            status = "OK"
                        else:
                            status = "FAILED"
                        out.write(f"{sample}: {status}\n")
                except FileNotFoundError:
                    out.write(f"{sample}: NO LOG\n")
