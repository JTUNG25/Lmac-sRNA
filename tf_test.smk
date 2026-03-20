#!/usr/bin/env python3

targetfinder = "docker pull quay.io/biocontainers/targetfinder:1.7--3"

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

        echo "Database: {input.db}" >> {log}
        echo "Score cutoff: {params.score}" >> {log}
        echo "" >> {log}

        # Check if database exists
        if [[ ! -f "{input.db}" ]]; then
            echo "ERROR: Database file not found: {input.db}" >> {log}
            echo "Available files in data/genome/:" >> {log}
            ls -la data/genome/ >> {log} 2>&1
            exit 1
        fi

        # Initialize output files
        mkdir -p results/tf
        for FMT in table gff json; do
            > results/tf/{wildcards.sample}_targets.$FMT
        done

        # Process sequences
        SRNA_NAME=""
        SRNA_SEQ=""
        PROCESSED_COUNT=0

        while IFS= read -r line || [[ -n "$line" ]]; do
            if [[ "$line" == ">"* ]]; then
                if [[ -n "$SRNA_NAME" && -n "$SRNA_SEQ" ]]; then
                    echo "Processing: $SRNA_NAME (length: ${{#SRNA_SEQ}})" >> {log}
                    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
                    
                    for FMT in table gff json; do
                        # CORRECTED COMMAND: Remove unsupported -f and -u flags
                        $TARGETFINDER_CMD -s "$SRNA_SEQ" -d {input.db} -q "$SRNA_NAME" -c {params.score} -p $FMT -t 1 -r >> results/tf/{wildcards.sample}_targets.$FMT 2>> {log}
                        
                        RETVAL=$?
                        if [ $RETVAL -ne 0 ]; then
                            echo "  WARNING: TargetFinder returned exit code $RETVAL for $SRNA_NAME ($FMT)" >> {log}
                        fi
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
            echo "Processing: $SRNA_NAME (length: ${{#SRNA_SEQ}})" >> {log}
            PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
            
            for FMT in table gff json; do
                $TARGETFINDER_CMD -s "$SRNA_SEQ" -d {input.db} -q "$SRNA_NAME" -c {params.score} -p $FMT -t 1 -r >> results/tf/{wildcards.sample}_targets.$FMT 2>> {log}
                
                RETVAL=$?
                if [ $RETVAL -ne 0 ]; then
                    echo "  WARNING: TargetFinder returned exit code $RETVAL for $SRNA_NAME ($FMT)" >> {log}
                fi
            done
        fi

        echo "" >> {log}
        echo "SUMMARY:" >> {log}
        echo "Processed $PROCESSED_COUNT sequences" >> {log}
        
        # Count results
        for FMT in table gff json; do
            HITS=$(grep -v '^#' results/tf/{wildcards.sample}_targets.$FMT 2>/dev/null | wc -l || echo 0)
            echo "$FMT hits: $HITS" >> {log}
        done
        
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
            out.write("=" * 50 + "\n\n")

            # Input stats
            with open(input.srna_fa) as f:
                srna_count = sum(1 for l in f if l.startswith(">"))

            # Output stats
            with open(input.table) as f:
                data_lines = [l for l in f if not l.startswith("#") and l.strip()]

            out.write(f"Input sRNA sequences: {srna_count}\n")
            out.write(f"Target hits found: {len(data_lines)}\n")

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
                    out.write(f"Mean score: {sum(scores)/len(scores):.2f}\n")
                    
                    out.write("\nScore distribution:\n")
                    for threshold in [1.0, 2.0, 3.0, 4.0, 5.0]:
                        count = sum(1 for s in scores if s <= threshold)
                        out.write(f"  <= {threshold}: {count} hits\n")
                        
                    out.write("\nTop 5 hits:\n")
                    sorted_hits = sorted(zip(scores, data_lines))[:5]
                    for i, (score, line) in enumerate(sorted_hits, 1):
                        parts = line.split("\t")
                        target = parts[1] if len(parts) > 1 else "Unknown"
                        out.write(f"  {i}. {parts[0]} -> {target[:50]}... (score: {score})\n")
            else:
                out.write("\nNO HITS FOUND\n")
                out.write("Possible reasons:\n")
                out.write("- Score cutoff too strict (try 5.0-6.0)\n")
                out.write("- Wrong database (need transcripts vs genome)\n")
                out.write("- sRNAs genuinely have no strong targets\n")
                out.write("- Database format issues\n")

rule summarise_hits:
    input:
        tables=expand("results/tf/{sample}_targets.table", sample=TEST_SAMPLES),
        diags=expand("results/tf/{sample}_diagnostics.txt", sample=TEST_SAMPLES),
    output:
        summary="results/tf/hit_summary.txt",
    run:
        total_hits = 0
        sample_results = {}
        
        with open(output.summary, "w") as out:
            out.write("TargetFinder Results Summary\n")
            out.write("=" * 40 + "\n")
            out.write(f"Database: {TARGET_DB}\n")
            out.write(f"Score cutoff: {SCORE_CUTOFF}\n")
            out.write(f"Samples analyzed: {len(TEST_SAMPLES)}\n\n")
            
            for table_file in input.tables:
                sample_name = table_file.split("/")[-1].replace("_targets.table", "")
                
                with open(table_file) as f:
                    data_lines = [l for l in f if not l.startswith("#") and l.strip()]
                
                hits = len(data_lines)
                total_hits += hits
                sample_results[sample_name] = hits
                
                # Get scores
                scores = []
                for line in data_lines:
                    parts = line.split("\t")
                    try:
                        scores.append(float(parts[4]))
                    except (IndexError, ValueError):
                        pass
                
                if scores:
                    out.write(f"{sample_name}: {hits} hits (scores: {min(scores):.1f}-{max(scores):.1f})\n")
                else:
                    out.write(f"{sample_name}: {hits} hits\n")
            
            out.write(f"\nTotal hits across all samples: {total_hits}\n")
            out.write(f"Average hits per sample: {total_hits/len(TEST_SAMPLES):.1f}\n")
            
            # Summary by hit count
            no_hits = sum(1 for hits in sample_results.values() if hits == 0)
            few_hits = sum(1 for hits in sample_results.values() if 1 <= hits <= 10)
            many_hits = sum(1 for hits in sample_results.values() if hits > 10)
            
            out.write(f"\nSample distribution:\n")
            out.write(f"  No hits: {no_hits} samples\n")
            out.write(f"  1-10 hits: {few_hits} samples\n")
            out.write(f"  >10 hits: {many_hits} samples\n")

rule workflow_diagnostics:
    input:
        logs=expand("logs/targetfinder/{sample}.log", sample=TEST_SAMPLES),
        tables=expand("results/tf/{sample}_targets.table", sample=TEST_SAMPLES),
    output:
        wf_diag="results/tf/workflow_diagnostics.txt",
    run:
        with open(output.wf_diag, "w") as out:
            out.write("TargetFinder Workflow Diagnostics\n")
            out.write("=" * 40 + "\n\n")
            out.write(f"Container: {targetfinder}\n")
            out.write(f"Database: {TARGET_DB}\n")
            out.write(f"Score cutoff: {SCORE_CUTOFF}\n\n")
            
            success_count = 0
            for log_file, sample in zip(input.logs, TEST_SAMPLES):
                try:
                    with open(log_file) as f:
                        content = f.read()
                        if "Found TargetFinder executable" in content and "Analysis complete" in content:
                            status = "SUCCESS"
                            success_count += 1
                        elif "ERROR" in content:
                            status = "ERROR"
                        else:
                            status = "INCOMPLETE"
                        out.write(f"{sample}: {status}\n")
                except FileNotFoundError:
                    out.write(f"{sample}: NO LOG FILE\n")
            
            out.write(f"\nSuccessful runs: {success_count}/{len(TEST_SAMPLES)}\n")
            
            if success_count == len(TEST_SAMPLES):
                out.write("\n✓ All samples processed successfully!\n")
                out.write("Check hit_summary.txt for results overview.\n")
            else:
                out.write(f"\n⚠ {len(TEST_SAMPLES) - success_count} samples failed.\n")
                out.write("Check individual log files for error details.\n")
