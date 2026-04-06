import gzip
from collections import defaultdict
from pathlib import Path
import csv

sample_files = list(Path("data/merged_srna").glob("*.merged.fastq.gz"))

length_rows = []
nuc_rows = []  # first nucleotide by length
pos_rows = []  # positional nucleotide (all positions, for lengths 20-24)

LENGTHS_OF_INTEREST = set(range(20, 27))  # adjust to your biology
NUCLEOTIDES = ["A", "T", "C", "G"]  # T in DNA-encoded fastq; swap to U in R

for f in sorted(sample_files):
    sample = f.name.replace(".merged.fastq.gz", "")

    length_counts = defaultdict(int)
    # first_nuc[length][nuc] = count
    first_nuc = defaultdict(lambda: defaultdict(int))
    # pos_nuc[length][position][nuc] = count
    pos_nuc = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with gzip.open(f, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seq = line.strip()
                l = len(seq)
                if 16 <= l <= 30:
                    length_counts[l] += 1
                    # First nucleotide
                    first_nuc[l][seq[0]] += 1
                    # Positional nucleotide (only for lengths of interest, capped at 24 nt)
                    if l in LENGTHS_OF_INTEREST:
                        for pos, nuc in enumerate(seq, start=1):
                            pos_nuc[l][pos][nuc] += 1

    # --- length counts ---
    for length, count in length_counts.items():
        length_rows.append({"sample": sample, "length": length, "count": count})

    # --- first nucleotide counts ---
    for length, nuc_counts in first_nuc.items():
        total = sum(nuc_counts.values())
        for nuc in NUCLEOTIDES:
            count = nuc_counts.get(nuc, 0)
            nuc_rows.append(
                {
                    "sample": sample,
                    "length": length,
                    "nucleotide": nuc,
                    "count": count,
                    "proportion": count / total if total > 0 else 0,
                }
            )

    # --- positional nucleotide counts ---
    for length, positions in pos_nuc.items():
        for pos, nuc_counts in positions.items():
            total = sum(nuc_counts.values())
            for nuc in NUCLEOTIDES:
                count = nuc_counts.get(nuc, 0)
                pos_rows.append(
                    {
                        "sample": sample,
                        "length": length,
                        "position": pos,
                        "nucleotide": nuc,
                        "count": count,
                        "proportion": count / total if total > 0 else 0,
                    }
                )

# --- write CSVs ---
with open("data/merged_srna/srna_length_counts.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["sample", "length", "count"])
    writer.writeheader()
    writer.writerows(length_rows)

with open("data/merged_srna/srna_first_nuc.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=["sample", "length", "nucleotide", "count", "proportion"]
    )
    writer.writeheader()
    writer.writerows(nuc_rows)

with open("data/merged_srna/srna_pos_nuc.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "sample",
            "length",
            "position",
            "nucleotide",
            "count",
            "proportion",
        ],
    )
    writer.writeheader()
    writer.writerows(pos_rows)

print("Done — wrote length, first-nucleotide, and positional nucleotide CSVs")
