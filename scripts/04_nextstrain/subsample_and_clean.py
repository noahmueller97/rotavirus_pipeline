"""
- Subsamples and cleans rotavirus datasets.
- Reduces large datasets to a manageable size while preserving temporal diversity and cleaning metadata fields.
- Used before phylogenetic reconstruction to ensure datasets remain computationally manageable.

Inputs:
- filtered sequence datasets
- metadata files

Outputs:
- subsampled FASTA files
- cleaned metadata tables
"""

from __future__ import annotations

import argparse
import csv
import random
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Set

from Bio import SeqIO


YEAR_RE = re.compile(r"(\d{4})")
DATE_YMD_RE = re.compile(r"^\s*(\d{4})(?:[-/](\d{1,2})(?:[-/](\d{1,2}))?)?\s*$")


def die(msg: str) -> None:
    raise SystemExit(f"[ERROR] {msg}")


def normalize_date(raw: str) -> str:
    """
    Return a Nextstrain-ish date string WITHOUT inventing precision:
      - YYYY-MM-DD -> keep (zero-pad)
      - YYYY-MM    -> keep (zero-pad)
      - YYYY       -> keep
      - otherwise  -> "" (or extract first year as fallback)
    """
    if raw is None:
        return ""
    s = str(raw).strip()
    if not s:
        return ""

    m = DATE_YMD_RE.match(s)
    if m:
        y, mo, d = m.group(1), m.group(2), m.group(3)
        if mo and d:
            return f"{y}-{int(mo):02d}-{int(d):02d}"
        if mo:
            return f"{y}-{int(mo):02d}"
        return y

    m2 = YEAR_RE.search(s)
    return m2.group(1) if m2 else ""


def year_from_date(norm_date: str) -> str:
    return norm_date[:4] if norm_date and len(norm_date) >= 4 else ""


def read_metadata(meta_tsv: Path) -> Tuple[List[str], List[Dict[str, str]]]:
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


def count_groups(rows: List[Dict[str, str]], col: str) -> Counter:
    c = Counter()
    for row in rows:
        v = (row.get(col) or "").strip()
        c[v if v else "MISSING"] += 1
    return c


def print_top(counter: Counter, n: int = 15, title: str = "") -> None:
    if title:
        print(title)
    for k, v in counter.most_common(n):
        print(f"  {k}\t{v}")
    if len(counter) > n:
        print(f"  ... ({len(counter) - n} more groups)")


def subsample_ids_total(groups: Dict[str, List[str]], total_n: int, rng: random.Random) -> Set[str]:
    """
    Distribute TOTAL total_n across groups proportional to group sizes.
    Deterministic given seed.
    """
    groups = {k: v for k, v in groups.items() if k and k != "MISSING" and v}
    if not groups:
        return set()

    sizes = {k: len(v) for k, v in groups.items()}
    total = sum(sizes.values())

    if total_n >= total:
        return set(x for v in groups.values() for x in v)

    alloc = {k: int(round(total_n * sizes[k] / total)) for k in groups}
    keys = sorted(groups.keys())

    
    for k in alloc:
        alloc[k] = max(0, min(alloc[k], sizes[k]))

    current = sum(alloc.values())

    
    while current < total_n:
        candidates = [(sizes[k] - alloc[k], k) for k in keys if alloc[k] < sizes[k]]
        if not candidates:
            break
        _, k = max(candidates)
        alloc[k] += 1
        current += 1

    while current > total_n:
        candidates = [(alloc[k], k) for k in keys if alloc[k] > 0]
        if not candidates:
            break
        _, k = max(candidates)
        alloc[k] -= 1
        current -= 1

    kept: Set[str] = set()
    for k in keys:
        ids = groups[k][:]
        rng.shuffle(ids)
        kept.update(ids[:alloc[k]])
    return kept


def subsample_ids_per_year(groups: Dict[str, List[str]], per_n: int, rng: random.Random) -> Set[str]:
    """
    Keep up to per_n sequences PER YEAR group.
    """
    groups = {k: v for k, v in groups.items() if k and k != "MISSING" and v}
    kept: Set[str] = set()
    for k in sorted(groups.keys()):
        ids = groups[k][:]
        rng.shuffle(ids)
        kept.update(ids[:min(per_n, len(ids))])
    return kept


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Normalize metadata date/year_group and subsample sequences deterministically."
    )
    ap.add_argument("--fasta-in", required=True, type=Path)
    ap.add_argument("--meta-in", required=True, type=Path)
    ap.add_argument("--fasta-out", required=True, type=Path)
    ap.add_argument("--meta-out", required=True, type=Path)

    ap.add_argument("--subsample", type=int, default=1000, help="Max sequences (TOTAL unless --subsample-per-year).")
    ap.add_argument("--subsample-per-year", action="store_true", help="If set: keep up to --subsample per year.")
    ap.add_argument("--seed", type=int, default=1)

    ap.add_argument("--year-col", default="year_group", help="Column name to store year grouping.")
    ap.add_argument("--drop-missing-year", action="store_true", help="Drop rows with no parseable year/date.")
    ap.add_argument("--strain-col", default="strain", help="Metadata column matching FASTA record IDs.")
    ap.add_argument("--date-col", default="date", help="Metadata date column (YYYY or YYYY-MM or YYYY-MM-DD).")

    args = ap.parse_args()

    if not args.fasta_in.exists():
        die(f"Missing FASTA: {args.fasta_in}")
    if not args.meta_in.exists():
        die(f"Missing metadata TSV: {args.meta_in}")

    
    fieldnames, rows = read_metadata(args.meta_in)
    if args.strain_col not in fieldnames:
        die(f"{args.meta_in} missing required column '{args.strain_col}'. Found: {fieldnames}")
    if args.date_col not in fieldnames:
        die(f"{args.meta_in} missing required column '{args.date_col}'. Found: {fieldnames}")

    out_fields = list(fieldnames)
    if args.year_col not in out_fields:
        out_fields.append(args.year_col)

    
    norm_rows: List[Dict[str, str]] = []
    dropped = 0
    for row in rows:
        raw_date = row.get(args.date_col, "")
        norm_date = normalize_date(raw_date)
        row[args.date_col] = norm_date

        y = year_from_date(norm_date)
        row[args.year_col] = y

        if args.drop_missing_year and not y:
            dropped += 1
            continue

        norm_rows.append(row)

    if not norm_rows:
        die("Metadata empty after date/year processing. Consider not using --drop-missing-year.")

    pre = count_groups(norm_rows, args.year_col)
    print_top(pre, n=15, title="Year groups BEFORE subsampling (top 15):")
    if args.drop_missing_year:
        print(f"Dropped missing-year rows: {dropped}")

    
    seq_dict = SeqIO.to_dict(SeqIO.parse(str(args.fasta_in), "fasta"))
    if not seq_dict:
        die(f"No sequences found in {args.fasta_in}")

    
    groups: Dict[str, List[str]] = defaultdict(list)
    for row in norm_rows:
        y = (row.get(args.year_col) or "").strip()
        sid = (row.get(args.strain_col) or "").strip()
        if not y or not sid:
            continue
        if sid in seq_dict: 
            groups[y].append(sid)

    rng = random.Random(args.seed)
    if args.subsample_per_year:
        kept_ids = subsample_ids_per_year(groups, args.subsample, rng)
        print(f"Subsampling mode: PER-YEAR up to {args.subsample}")
    else:
        kept_ids = subsample_ids_total(groups, args.subsample, rng)
        print(f"Subsampling mode: TOTAL {args.subsample} distributed across years")

    if not kept_ids:
        die("Subsampling selected 0 sequences. Check year parsing, strain IDs, or --drop-missing-year.")

    
    args.fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with args.fasta_out.open("w") as f:
        for sid in sorted(kept_ids):
            rec = seq_dict.get(sid)
            if rec is not None:
                SeqIO.write(rec, f, "fasta")

    
    sub_rows = [row for row in norm_rows if (row.get(args.strain_col) or "").strip() in kept_ids]
    if not sub_rows:
        die("Subsampled metadata ended up empty (strain IDs mismatch FASTA headers).")

    write_metadata(args.meta_out, out_fields, sub_rows)

    post = count_groups(sub_rows, args.year_col)
    print_top(post, n=15, title="Year groups AFTER subsampling (top 15):")
    print(f"Wrote FASTA: {args.fasta_out}  (n={len(kept_ids)})")
    print(f"Wrote META : {args.meta_out}  (n={len(sub_rows)})")


if __name__ == "__main__":
    main()