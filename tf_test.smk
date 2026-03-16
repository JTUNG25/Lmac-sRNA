#!/usr/bin/env python3

targetfinder = "docker://quay.io/biocontainers/targetfinder:1.7--0"

SCORE_CUTOFF = 4.0
TARGET_DB    = "data/genome/JN3_transcript_clean.fa"
TEST_SAMPLES = ["D1_top50", "A3_top50", "R1_top50"]

rule target:
    """Final rule - generate all outputs including diagnostics"""
    input:
        expand(
            "results/test/{sample}_targets.{fmt}",
            sample=TEST_SAMPLES,
            fmt=["table", "gff", "json"],
        ),
        expand(
            "results/test/{sample}_diagnostics.txt",
            sample=TEST_SAMPLES,
        ),
        "results/test/hit_summary.txt",
        "results/test/workflow_diagnostics.txt",

rule targetfinder:
    """Run TargetFinder on each sRNA sample with multiple output formats"""
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
    log:
        "logs/targetfinder/{sample}.log",
    shell:
        """
        echo "========================================" > {log}
        echo "TargetFinder Analysis: {wildcards.sample}" >> {log}
        echo "========================================" >> {log}
        echo "Input FASTA: {input.srna_fa}" >> {log}
        echo "Database: {input.db}" >> {log}
        echo "Score cutoff: {SCORE_CUTOFF}" >> {log}
        echo "Start time: $(date)" >> {log}
        echo "" >> {log}
        
        # Diagnostic: Check input file
        echo "Input file check:" >> {log}
        wc -l {input.srna_fa} >> {log}
        head -4 {input.srna_fa} >> {log}
        echo "" >> {log}
        
        # Count sequences
        SRNA_COUNT=$(grep -c "^>" {input.srna_fa})
        echo "Number of sRNA sequences: $SRNA_COUNT" >> {log}
        echo "" >> {log}
        
        # Extract sequences for TargetFinder
        SRNA_SEQ=$(grep -v '^>' {input.srna_fa} | tr -d '\n')
        SRNA_NAME="{wildcards.sample}"
        
        echo "Processing $SRNA_COUNT sequences..." >> {log}
        echo "" >> {log}
        
        # Run TargetFinder for each format
        for FMT in table gff json; do
            echo "Running TargetFinder format: $FMT" >> {log}
            
            targetfinder.pl \
                -s "$SRNA_SEQ" \
                -d {input.db} \
                -q "$SRNA_NAME" \
                -c {SCORE_CUTOFF} \
                -p $FMT \
                > results/test/{wildcards.sample}_targets.$FMT \
                2>> {log}
            
            echo "Output lines for $FMT: $(wc -l < results/test/{wildcards.sample}_targets.$FMT)" >> {log}
        done
        
        echo "" >> {log}
        echo "End time: $(date)" >> {log}
        echo "Status: SUCCESS" >> {log}
        """

rule diagnose_sample:
    """Generate detailed diagnostic report for each sample"""
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
            
            # Input statistics
            with open(input.srna_fa) as f:
                srna_lines = f.readlines()
                srna_count = sum(1 for l in srna_lines if l.startswith(">"))
                seq_count = sum(1 for l in srna_lines if not l.startswith(">") and l.strip())
            
            out.write("INPUT STATISTICS:\n")
            out.write("-" * 60 + "\n")
            out.write(f"Input file: {input.srna_fa}\n")
            out.write(f"sRNA sequences: {srna_count}\n")
            out.write(f"Sequence lines: {seq_count}\n")
            out.write(f"Total file lines: {len(srna_lines)}\n\n")
            
            # Output analysis
            with open(input.table) as f:
                table_lines = f.readlines()
            
            comment_lines = sum(1 for l in table_lines if l.startswith("#"))
            data_lines = sum(1 for l in table_lines if not l.startswith("#") and l.strip())
            
            out.write("OUTPUT STATISTICS:\n")
            out.write("-" * 60 + "\n")
            out.write(f"Total table lines: {len(table_lines)}\n")
            out.write(f"Comment lines: {comment_lines}\n")
            out.write(f"Data lines (predictions): {data_lines}\n\n")
            
            # Parse scores
            scores = []
            targets = []
            for line in table_lines:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                try:
                    score = float(parts[4])
                    scores.append(score)
                    targets.append(parts)
                except (IndexError, ValueError):
                    pass
            
            out.write("PREDICTION STATISTICS:\n")
            out.write("-" * 60 + "\n")
            
            if scores:
                out.write(f"Total predictions: {len(scores)}\n")
                out.write(f"Score range: {min(scores):.2f} - {max(scores):.2f}\n")
                out.write(f"Mean score: {sum(scores)/len(scores):.2f}\n")
                out.write(f"Median score: {sorted(scores)[len(scores)//2]:.2f}\n\n")
                
                out.write("SCORE DISTRIBUTION:\n")
                for threshold in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]:
                    count = sum(1 for s in scores if s <= threshold)
                    pct = (count / len(scores)) * 100
                    out.write(f"  <= {threshold:.1f}: {count:4d} predictions ({pct:5.1f}%)\n")
                
                out.write("\nTOP 5 PREDICTIONS:\n")
                for i, (score, target) in enumerate(sorted(zip(scores, targets), key=lambda x: x[0])[:5]):
                    out.write(f"  {i+1}. Score: {score:.2f} - Target: {target[0] if target else 'N/A'}\n")
            else:
                out.write("NO PREDICTIONS FOUND\n")
                out.write("Interpretation:\n")
                out.write("  - Score cutoff may be too strict\n")
                out.write("  - sRNAs may not have good matches in database\n")
                out.write("  - Try with cutoff=5.0 or higher\n")
            
            out.write("\n" + "=" * 60 + "\n")
            out.write("INTERPRETATION GUIDE:\n")
            out.write("=" * 60 + "\n")
            out.write("0 hits      : Cutoff too strict, try 5.0 or 6.0\n")
            out.write("1-10 hits   : Excellent! High-confidence targets\n")
            out.write("10-50 hits  : Good, selective targeting\n")
            out.write("50-100 hits : Reasonable, moderate specificity\n")
            out.write(">100 hits   : May want tighter cutoff\n")
            out.write("\n")

rule summarise_hits:
    """Summarize all hits across samples"""
    input:
        tables = expand("results/test/{sample}_targets.table", sample=TEST_SAMPLES),
        diags = expand("results/test/{sample}_diagnostics.txt", sample=TEST_SAMPLES),
    output:
        summary = "results/test/hit_summary.txt",
    run:
        with open(output.summary, "w") as out:
            out.write("=" * 60 + "\n")
            out.write("TargetFinder Test Run - Overall Summary\n")
            out.write("=" * 60 + "\n\n")
            
            out.write("PARAMETERS:\n")
            out.write("-" * 60 + "\n")
            out.write(f"Score cutoff  : {SCORE_CUTOFF}\n")
            out.write(f"Database      : {TARGET_DB}\n")
            out.write(f"Samples tested: {', '.join(TEST_SAMPLES)}\n\n")
            
            out.write("RESULTS PER SAMPLE:\n")
            out.write("-" * 60 + "\n")
            
            total_hits = 0
            all_scores = []
            
            for table_file in input.tables:
                sample_name = table_file.split("/")[-1].replace("_targets.table", "")
                
                with open(table_file) as f:
                    lines = [l for l in f if not l.startswith("#") and l.strip()]
                
                hits = len(lines)
                total_hits += hits
                
                scores = []
                for line in lines:
                    parts = line.split("\t")
                    try:
                        scores.append(float(parts[4]))
                        all_scores.append(float(parts[4]))
                    except (IndexError, ValueError):
                        pass
                
                out.write(f"\n{sample_name}:\n")
                out.write(f"  Hits: {hits}\n")
                if scores:
                    out.write(f"  Score range: {min(scores):.2f} - {max(scores):.2f}\n")
                    out.write(f"  Mean score: {sum(scores)/len(scores):.2f}\n")
                else:
                    out.write(f"  Status: No targets found\n")
            
            out.write(f"\n\nOVERALL STATISTICS:\n")
            out.write("-" * 60 + "\n")
            out.write(f"Total samples: {len(TEST_SAMPLES)}\n")
            out.write(f"Total hits: {total_hits}\n")
            
            if all_scores:
                out.write(f"Overall score range: {min(all_scores):.2f} - {max(all_scores):.2f}\n")
                out.write(f"Overall mean score: {sum(all_scores)/len(all_scores):.2f}\n")
            
            out.write("\n" + "=" * 60 + "\n")
            out.write("NEXT STEPS:\n")
            out.write("=" * 60 + "\n")
            
            if total_hits == 0:
                out.write("No targets found across all samples.\n")
                out.write("Options:\n")
                out.write("  1. Increase score cutoff (e.g., 5.0, 6.0)\n")
                out.write("  2. Check if sRNAs match genome database\n")
                out.write("  3. Verify sequence formats are correct\n")
            else:
                out.write(f"Found {total_hits} predictions across samples.\n")
                out.write("Next steps:\n")
                out.write("  1. Review individual sample diagnostics\n")
                out.write("  2. Check score distribution\n")
                out.write("  3. Decide if cutoff is appropriate\n")
                out.write("  4. Run full workflow on all samples\n")
            
            out.write("=" * 60 + "\n")

rule workflow_diagnostics:
    """Generate workflow-level diagnostic report"""
    input:
        tables = expand("results/test/{sample}_targets.table", sample=TEST_SAMPLES),
        logs = expand("logs/targetfinder/{sample}.log", sample=TEST_SAMPLES),
    output:
        wf_diag = "results/test/workflow_diagnostics.txt",
    run:
        with open(output.wf_diag, "w") as out:
            out.write("=" * 70 + "\n")
            out.write("WORKFLOW DIAGNOSTIC REPORT\n")
            out.write("=" * 70 + "\n\n")
            
            out.write("EXECUTION INFORMATION:\n")
            out.write("-" * 70 + "\n")
            out.write(f"Timestamp: {__import__('datetime').datetime.now()}\n")
            out.write(f"Container: {targetfinder}\n")
            out.write(f"Score cutoff: {SCORE_CUTOFF}\n")
            out.write(f"Database: {TARGET_DB}\n")
            out.write(f"Samples: {len(TEST_SAMPLES)}\n\n")
            
            out.write("LOG FILE SUMMARIES:\n")
            out.write("-" * 70 + "\n")
            for log_file, sample in zip(input.logs, TEST_SAMPLES):
                out.write(f"\n{sample} log:\n")
                try:
                    with open(log_file) as f:
                        log_content = f.read()
                        # Extract key info
                        if "SUCCESS" in log_content:
                            out.write("  Status: SUCCESS ✓\n")
                        elif "ERROR" in log_content:
                            out.write("  Status: ERROR ✗\n")
                        else:
                            out.write("  Status: UNKNOWN\n")
                        
                        # Count sequences processed
                        for line in log_content.split("\n"):
                            if "sRNA sequences:" in line:
                                out.write(f"  {line.strip()}\n")
                except FileNotFoundError:
                    out.write("  Log file not found\n")
            
            out.write("\n" + "=" * 70 + "\n")
            out.write("RECOMMENDATIONS:\n")
            out.write("=" * 70 + "\n")
            out.write("1. Check individual sample diagnostic files\n")
            out.write("2. Review hit_summary.txt for overall statistics\n")
            out.write("3. Adjust score cutoff if needed\n")
            out.write("4. If satisfied, run full workflow on grouped samples\n")
            out.write("=" * 70 + "\n")