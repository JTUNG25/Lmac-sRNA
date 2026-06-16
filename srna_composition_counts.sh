#!/bin/bash
#SBATCH --job-name=srna_comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/srna_comp_%j.log

source ~/.bashrc
conda activate srna_tools

cd /QRISdata/Q9140/lmac/lmac_srna

export TMPDIR=/QRISdata/Q9140/lmac/lmac_srna/tmp
mkdir -p "$TMPDIR"

INDIR="data/fastp_srna"
OUTDIR="data/fastp_srna"

# write headers
echo "sample,length,count"                               > "$OUTDIR/srna_length_counts.csv"
echo "sample,length,nucleotide,count,proportion"         > "$OUTDIR/srna_first_nuc.csv"
echo "sample,length,position,nucleotide,count,proportion" > "$OUTDIR/srna_pos_nuc.csv"

for f in "$INDIR"/*.merged.fastq.gz; do
    sample=$(basename "$f" .merged.fastq.gz)
    echo "Processing $sample..."

    # --- length counts + first nucleotide in one pass ---
    seqkit seq -m 16 -M 30 -j 8 "$f" 2>/dev/null | \
    seqkit fx2tab -l -Q --only-id 2>/dev/null | \
    awk -v s="$sample" \
        -v len_out="$OUTDIR/srna_length_counts.csv" \
        -v nuc_out="$OUTDIR/srna_first_nuc.csv" '
    {
        seq = $2
        len = int($3)
        if (len < 16 || len > 30) next

        # length counts
        length_count[len]++

        # first nucleotide
        first = substr(seq, 1, 1)
        if (first == "A" || first == "T" || first == "C" || first == "G") {
            first_count[len][first]++
            first_total[len]++
        }
    }
    END {
        for (len in length_count) {
            print s "," len "," length_count[len] >> len_out
        }
        for (len in first_count) {
            for (nuc in first_count[len]) {
                count = first_count[len][nuc]
                prop  = count / first_total[len]
                print s "," len "," nuc "," count "," prop >> nuc_out
            }
        }
    }'

    # --- positional nucleotide ---
    seqkit seq -m 18 -M 21 -j 8 "$f" 2>/dev/null | \
    seqkit fx2tab -l -Q --only-id 2>/dev/null | \
    awk -v s="$sample" \
        -v pos_out="$OUTDIR/srna_pos_nuc.csv" '
    {
        seq = $2
        len = int($3)
        for (i = 1; i <= len; i++) {
            nuc = substr(seq, i, 1)
            if (nuc == "A" || nuc == "T" || nuc == "C" || nuc == "G") {
                pos_count[len][i][nuc]++
                pos_total[len][i]++
            }
        }
    }
    END {
        for (len in pos_count) {
            for (pos in pos_count[len]) {
                for (nuc in pos_count[len][pos]) {
                    count = pos_count[len][pos][nuc]
                    prop  = count / pos_total[len][pos]
                    print s "," len "," pos "," nuc "," count "," prop >> pos_out
                }
            }
        }
    }'

    echo "  Done: $sample"
done

echo "All samples complete"
