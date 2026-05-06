#!/usr/bin/env python3

# =============================================================================
# MiRNATarget (psRNATarget reimplementation) Snakemake Pipeline
# Run from: /QRISdata/Q9140/lmac/lmac_srna/
# Usage:    snakemake -s miRNATarget.smk --profile profiles/bunya/
# =============================================================================

FASTA36_DIR     = "/QRISdata/Q9140/lmac/lmac_srna/fasta36"
SSEARCH36       = f"{FASTA36_DIR}/bin/ssearch36"
MIRNATARGET_DIR = "/QRISdata/Q9140/lmac/lmac_srna/MiRNATarget"


# --- Parameters ---------------------------------------------------------------
E_CUTOFF        = 4.0   # Max penalty score — 4.0 matches TargetFinder for comparison
                         # psRNATarget default is 5.0
SEED_START      = 2     # Default psRNATarget seed start
SEED_END        = 13    # Default psRNATarget seed end (TargetFinder uses 12)
PENALTY_MULT    = 1.5   # Multiplier applied to mismatches in seed region only
MAX_SEED_MM     = 2     # Max mismatches allowed in seed region
GAP_CUTOFF      = 1     # Max gaps allowed in alignment
TOTAL_MM_CUTOFF = 8     # Max total mismatches allowed
GU_CUTOFF       = 7     # Max G:U wobble pairs allowed
HSP_CUTOFF      = 14    # Min alignment length (HSP) to keep
MAX_ALIGN_LEN   = 22    # Max alignment length to keep
ALIGN_LENGTH    = 19    # How many bp from 5' of sRNA are scored

TARGET_DB = "data/genome/JN3_transcript_clean.fa"

TEST_SAMPLES = [
    "A1-1_top50",
    "A1-2_top50",
    "A1-3_top50",
    "A3_top50",
    "A13-1_top50",
    "A13-2_top50",
    "D1_top50",
    "D2-2_top50",
    "D2-3_top50",
    "D2-4_top50",
    "R1_top50",
    "R2-2_top50",
    "R2-3_top50",
    "R2-4_top50",
    "R2-5_top50",
    "R12-1_top50",
    "R12-2_top50",
    "R12-3_top50",
]


# =============================================================================
# Target rule
# =============================================================================
rule target:
    input:
        expand(
            "results/mirnatarget/{sample}_targets.tsv",
            sample=TEST_SAMPLES,
        ),
        "results/mirnatarget/hit_summary.txt",


# =============================================================================
# Step 1: Reverse complement the sRNA fasta
# SSEARCH aligns the revcomp of the sRNA against the forward transcript strand,
# reproducing psRNATarget's internal strand handling.
# =============================================================================

rule revcomp_srna:
    input:
        srna_fa = "data/all/edger_seq/{sample}.fasta",
    output:
        rc_fa = temp("results/mirnatarget/{sample}_rc.fasta"),
    shell:
        """
        seqkit seq --reverse --complement --rna2dna {input.srna_fa} \
            > {output.rc_fa}
        """

# =============================================================================
# Step 2: Run SSEARCH36 alignment
# Parameters validated by MiRNATarget authors to reproduce psRNATarget output.
# -f -8 -g -3  : gap open / extend penalties
# -E 10000     : keep up to 10000 hits per query
# -b 200       : max alignments reported per query
# -r +4/-3     : nucleotide match reward / mismatch penalty
# -n           : treat as nucleotide
# -U           : RNA alphabet (U not T)
# -W 10        : word size
# -N 20000     : max DB sequence length per chunk
# =============================================================================

rule ssearch_align:
    input:
        rc_fa = "results/mirnatarget/{sample}_rc.fasta",
        db    = TARGET_DB,
    output:
        ssearch_out = temp("results/mirnatarget/{sample}_ssearch.txt"),
    threads: 8
    resources:
        mem_mb  = 16000,
        runtime = 120,
    shell:
        """
        {SSEARCH36} \
            -f -8 -g -3 \
            -E 10000 \
            -T {threads} \
            -b 200 \
            -r +4/-3 \
            -n -U \
            -W 10 \
            -N 20000 \
            {input.rc_fa} \
            {input.db} \
            > {output.ssearch_out}
        """


# =============================================================================
# Step 3: Parse and score alignments using psRNATarget scheme
# Pipes parse_ssearch.py -> parse_mirna_targets.py to avoid intermediate files.
# =============================================================================
rule mirnatarget_parse:
    input:
        ssearch_out = "results/mirnatarget/{sample}_ssearch.txt",
    output:
        tsv = "results/mirnatarget/{sample}_targets.tsv",
    params:
        parse_ssearch   = f"{MIRNATARGET_DIR}/parse_ssearch.py",
        parse_targets   = f"{MIRNATARGET_DIR}/parse_mirna_targets.py",
        e_cutoff        = E_CUTOFF,
        seed_start      = SEED_START,
        seed_end        = SEED_END,
        penalty_mult    = PENALTY_MULT,
        max_seed_mm     = MAX_SEED_MM,
        gap_cutoff      = GAP_CUTOFF,
        total_mm_cutoff = TOTAL_MM_CUTOFF,
        gu_cutoff       = GU_CUTOFF,
        hsp_cutoff      = HSP_CUTOFF,
        max_align_len   = MAX_ALIGN_LEN,
        align_length    = ALIGN_LENGTH,
    shell:
        """
        python3 {params.parse_ssearch} -i {input.ssearch_out} \
            | python3 {params.parse_targets} \
                --E_cutoff                 {params.e_cutoff}        \
                --seed_start               {params.seed_start}      \
                --seed_end                 {params.seed_end}        \
                --penalty_multiplier       {params.penalty_mult}    \
                --num_mismatch_seed        {params.max_seed_mm}     \
                --gap_cutoff               {params.gap_cutoff}      \
                --total_mismatches_cutoff  {params.total_mm_cutoff} \
                --GUs_cutoff               {params.gu_cutoff}       \
                --hsp_cutoff               {params.hsp_cutoff}      \
                --maximum_alignment_length {params.max_align_len}   \
                --alignment_length         {params.align_length}    \
            > {output.tsv}
        """


# =============================================================================
# Step 4: Summarise hits across all samples
# Output columns from parse_mirna_targets.py (verified from source code):
#   0: query_id    1: subject_id   2: match_aln    3: query_aln
#   4: subject_aln 5: q_start      6: q_end        7: s_start
#   8: s_end       9: expect_value 10: strand
# Score (expect_value) is at index 9 — NOT index 2
# =============================================================================
rule summarise_hits:
    input:
        tsvs = expand("results/mirnatarget/{sample}_targets.tsv", sample=TEST_SAMPLES),
    output:
        summary = "results/mirnatarget/hit_summary.txt",
    run:
        sample_results = {}

        for tsv_file in input.tsvs:
            sample_name = tsv_file.split("/")[-1].replace("_targets.tsv", "")
            with open(tsv_file) as f:
                data_lines = [
                    l for l in f
                    if not l.startswith("#") and l.strip()
                ]
            scores = []
            for line in data_lines:
                parts = line.split("\t")
                try:
                    scores.append(float(parts[9]))   # expect_value = col index 9
                except (IndexError, ValueError):
                    pass
            sample_results[sample_name] = scores

        total_hits = sum(len(s) for s in sample_results.values())

        with open(output.summary, "w") as out:
            out.write("MiRNATarget (psRNATarget) Results Summary\n")
            out.write("=" * 45 + "\n")
            out.write(f"Database         : {TARGET_DB}\n")
            out.write(f"E_cutoff         : {E_CUTOFF}\n")
            out.write(f"Seed region      : pos {SEED_START}-{SEED_END}\n")
            out.write(f"Penalty mult     : {PENALTY_MULT}x (seed mismatches only)\n")
            out.write(f"Samples          : {len(TEST_SAMPLES)}\n\n")

            for sample_name, scores in sample_results.items():
                hits = len(scores)
                if scores:
                    out.write(
                        f"{sample_name}: {hits} hits  "
                        f"(scores {min(scores):.2f} - {max(scores):.2f})\n"
                    )
                else:
                    out.write(f"{sample_name}: 0 hits\n")

            out.write(f"\nTotal hits          : {total_hits}\n")
            out.write(
                f"Average per sample  : {total_hits / len(TEST_SAMPLES):.1f}\n"
            )

            no_hits   = sum(1 for s in sample_results.values() if len(s) == 0)
            few_hits  = sum(1 for s in sample_results.values() if 1 <= len(s) <= 10)
            many_hits = sum(1 for s in sample_results.values() if len(s) > 10)

            out.write(f"\nSample distribution:\n")
            out.write(f"  No hits  : {no_hits}\n")
            out.write(f"  1-10 hits: {few_hits}\n")
            out.write(f"  >10 hits : {many_hits}\n")

            out.write(f"""
Notes
-----
Database: transcript sequences (JN3_transcript_clean.fa), same as TargetFinder.
E_cutoff set to {E_CUTOFF} to match TargetFinder for cross-tool comparison (psRNATarget default: 5.0).
Seed region: positions {SEED_START}-{SEED_END} (TargetFinder uses 2-12).
Seed mismatches penalised with {PENALTY_MULT}x multiplier; G:U pairs (0.5) are NOT multiplied.
Gap opening: 2.0; gap extension: 0.5; extra query-gap penalty: 1.0.
Score column in output TSV: column 10 (0-based index 9), header 'expect_value'.
""")
