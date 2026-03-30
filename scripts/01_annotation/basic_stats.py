"""
- Computes basic statistics for annotated rotavirus sequences.
- Generates summary statistics for annotated sequences, primarily used to compare segment counts and dataset composition.
- Used for quick diagnostic statistics after the annotation step in the pipeline.

Inputs:
- annotated_genbank.ndjson

Outputs:
- segment_comparison.tsv

"""

import argparse
import json
from collections import Counter

def norm(x):
    if not x:
        return "unknown"
    return str(x).strip().lower()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--infile", default="results/annotated_genbank.ndjson")
    p.add_argument("--outfile", default="results/00_primary/segment_comparison.tsv")
    args = p.parse_args()

    INFILE = args.infile
    OUTFILE = args.outfile

    ncbi_counts = Counter()
    diamond_counts = Counter()
    agree_counts = Counter()

    with open(INFILE) as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)

            ncbi = norm(r.get("Segment"))
            diamond = norm(r.get("new_segment"))

            ncbi_counts[ncbi] += 1
            diamond_counts[diamond] += 1

            if ncbi == diamond and ncbi != "unknown":
                agree_counts[ncbi] += 1

    segments = sorted(set(ncbi_counts) | set(diamond_counts))

  
    import os
    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

    with open(OUTFILE, "w") as out:
        out.write("segment\tncbi_count\tdiamond_count\tagree_count\n")
        for seg in segments:
            out.write(
                f"{seg}\t"
                f"{ncbi_counts.get(seg,0)}\t"
                f"{diamond_counts.get(seg,0)}\t"
                f"{agree_counts.get(seg,0)}\t"
                f"{agree_counts.get(seg,0)}\n"
            )

    print(f"Wrote {OUTFILE}")

if __name__ == "__main__":
    main()