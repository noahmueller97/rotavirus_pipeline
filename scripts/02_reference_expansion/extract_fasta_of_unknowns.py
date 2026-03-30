"""
- Extracts FASTA sequences for unannotated rotavirus records.
- Reads the annotated NDJSON dataset and extracts sequences without confident segment annotation.
- These sequences are used for MMseqs clustering in order to identify potential missing reference segments.

Inputs:
- annotated_genbank.ndjson

Outputs:
- diamond_unknowns.fasta

"""

from __future__ import annotations
import argparse
import json
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", default="results/annotated_genbank.ndjson")
    ap.add_argument("--outfasta", default="results/01_unknowns/diamond_unknowns.fasta")
    args = ap.parse_args()

    infile = Path(args.infile)
    outfasta = Path(args.outfasta)
    outfasta.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with infile.open() as f, outfasta.open("w") as out:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)

            seg = (r.get("new_segment") or "").strip().lower()
            
            if seg and seg != "unknown":
                continue

            acc = (r.get("Accession") or "").strip()
            seq = (r.get("sequence") or "").strip().replace("\n", "")
            if not acc or not seq:
                continue

            out.write(f">{acc}\n{seq}\n")
            n += 1

    print(f"Wrote {n} sequences to {outfasta}")

if __name__ == "__main__":
    main()