"""
- Exports size statistics for MMseqs clusters.
- Reads clustering output from MMseqs and computes cluster size statistics, typically to identify large clusters of unknown sequences.
- Used to guide reference expansion by highlighting frequently occurring unknown sequence clusters.

Inputs:
- MMseqs cluster TSV file

Outputs:
- unknown_cluster_sizes_round2.tsv

"""
from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster-tsv", default="results/01_unknowns/mmseqs_unknowns_round2/unk2_cluster.tsv")
    ap.add_argument("--out", default="results/01_unknowns/unknown_cluster_sizes_round2.tsv")
    args = ap.parse_args()

    infile = Path(args.cluster_tsv)
    outfile = Path(args.out)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    counts = Counter()
    with infile.open() as f:
        for line in f:
            if not line.strip():
                continue
            rep = line.rstrip("\n").split("\t")[0]
            counts[rep] += 1

    with outfile.open("w") as out:
        out.write("representative\tcluster_size\n")
        for rep, n in counts.most_common():
            out.write(f"{rep}\t{n}\n")

    print(f"Wrote {len(counts)} clusters to {outfile}")

if __name__ == "__main__":
    main()