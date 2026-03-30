"""
- Filters alignments based on quality control criteria.
- Removes sequences from alignments that fail QC criteria, such as excessive gaps or poor coverage.

Inputs:
- aligned FASTA datasets

Outputs:
- filtered alignments suitable for phylogenetic analysis
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Dict, Set

from Bio import SeqIO


def die(msg: str) -> None:
    raise SystemExit(f"[ERROR] {msg}")


def read_metadata(meta_tsv: Path) -> tuple[List[str], List[Dict[str, str]]]:
    with meta_tsv.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            die(f"{meta_tsv} has no header row")
        rows = [row for row in r]
        return list(r.fieldnames), rows


def write_metadata(meta_tsv: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    meta_tsv.parent.mkdir(parents=True, exist_ok=True)
    with meta_tsv.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Filter an aligned FASTA by coverage/N fraction and filter metadata to matching strain IDs."
    )
    ap.add_argument("--aligned-fasta", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--out-aligned", required=True, type=Path)
    ap.add_argument("--out-meta", required=True, type=Path)
    ap.add_argument("--min-coverage", type=float, default=0.80, help="Min non-gap fraction (1 - gap_fraction).")
    ap.add_argument("--max-n-frac", type=float, default=0.05, help="Max N fraction among non-gap sites.")
    ap.add_argument("--strain-col", default="strain", help="Metadata column matching FASTA record IDs.")
    args = ap.parse_args()

    if not args.aligned_fasta.exists():
        die(f"Missing aligned FASTA: {args.aligned_fasta}")
    if not args.metadata.exists():
        die(f"Missing metadata TSV: {args.metadata}")

    recs = list(SeqIO.parse(str(args.aligned_fasta), "fasta"))
    if not recs:
        die(f"Empty alignment: {args.aligned_fasta}")

    L = len(str(recs[0].seq))
    if L == 0:
        die(f"Alignment has zero length: {args.aligned_fasta}")

    kept_recs = []
    kept_ids: Set[str] = set()

    for r in recs:
        s = str(r.seq).upper()
        gaps = s.count("-")
        nongap = L - gaps
        coverage = (nongap / L) if L else 0.0
        n_frac = (s.count("N") / nongap) if nongap else 1.0

        if coverage >= args.min_coverage and n_frac <= args.max_n_frac:
            kept_recs.append(r)
            kept_ids.add(r.id)

    if not kept_recs:
        die(
            f"QC removed everything (min_coverage={args.min_coverage}, max_n_frac={args.max_n_frac}). "
            f"Try relaxing thresholds."
        )

    args.out_aligned.parent.mkdir(parents=True, exist_ok=True)
    with args.out_aligned.open("w") as f:
        SeqIO.write(kept_recs, f, "fasta")

    fieldnames, rows = read_metadata(args.metadata)
    if args.strain_col not in fieldnames:
        die(f"{args.metadata} missing column '{args.strain_col}'. Found: {fieldnames}")

    out_rows = []
    for row in rows:
        sid = (row.get(args.strain_col) or "").strip()
        if sid in kept_ids:
            out_rows.append(row)

    if not out_rows:
        die("After QC, metadata became empty (strain IDs mismatch?).")

    write_metadata(args.out_meta, fieldnames, out_rows)

    print(f"[QC] kept_seqs={len(kept_recs)} / total={len(recs)}")
    print(f"[QC] wrote aligned: {args.out_aligned}")
    print(f"[QC] wrote meta   : {args.out_meta} (rows={len(out_rows)})")


if __name__ == "__main__":
    main()