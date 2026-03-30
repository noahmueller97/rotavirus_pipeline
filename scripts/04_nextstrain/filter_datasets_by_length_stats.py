"""
- Filters rotavirus datasets based on sequence length statistics.
- Uses previously computed segment length statistics to remove sequences that are substantially shorter than expected.

Inputs:
- sequence length statistics
- segment datasets

Outputs:
- filtered datasets for downstream analysis
"""

from __future__ import annotations

from pathlib import Path
import argparse
import csv
import math

def read_length_stats(stats_tsv: Path) -> dict:
    with stats_tsv.open("r") as f:
        header = f.readline().strip().split("\t")
        line = f.readline().strip().split("\t")
    if len(header) < 5 or len(line) < 5:
        raise ValueError(f"Bad stats file: {stats_tsv}")
    d = dict(zip(header, line))

    out = {}
    for k in ["n", "min", "p05", "p25", "median", "p75", "p95", "max", "mean"]:
        if k in d and d[k] != "":
            out[k] = float(d[k]) if k != "n" else int(float(d[k]))
    return out

def fasta_records(path: Path):
    seq_id = None
    chunks = []
    with path.open("r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].split()[0].split("|")[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if seq_id is not None:
            yield seq_id, "".join(chunks)

def filter_fasta(in_fasta: Path, out_fasta: Path, min_len: int) -> set[str]:
    kept_ids = set()
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w") as out:
        for sid, seq in fasta_records(in_fasta):
            seq2 = seq.replace("-", "")
            if len(seq2) >= min_len:
                kept_ids.add(sid)
                out.write(f">{sid}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + "\n")
    return kept_ids

def filter_metadata(in_tsv: Path, out_tsv: Path, kept_ids: set[str]) -> int:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    kept_rows = 0
    with in_tsv.open("r", newline="") as fin, out_tsv.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header in metadata: {in_tsv}")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        id_cols = [c for c in ["strain", "name", "seq_id", "accession", "Accession"] if c in reader.fieldnames]
        if not id_cols:
            id_cols = [reader.fieldnames[0]]

        for row in reader:
            row_id = None
            for c in id_cols:
                v = row.get(c, "")
                if v:
                    row_id = v.split()[0]
                    break
            if row_id and row_id in kept_ids:
                writer.writerow(row)
                kept_rows += 1
    return kept_rows

def write_filter_params(path: Path, params: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["key", "value"])
        for k, v in params.items():
            w.writerow([k, v])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets-root", type=Path, default=Path("results/02_datasets"))
    ap.add_argument("--nextstrain-root", type=Path, default=Path("results/03_nextstrain"))
    ap.add_argument("--median-frac", type=float, default=0.80,
                    help="min_len = max(p05, floor(median*median-frac))")
    ap.add_argument("--minlen-floor", type=int, default=0,
                    help="absolute floor for min_len (0 disables)")
    ap.add_argument("--only", type=str, default="",
                    help="optional filter: 'A/vp4' or 'A/*' or '*/vp4'")
    args = ap.parse_args()

    ds_root = args.datasets_root
    filtered_root = args.nextstrain_root / "01_filtered"

    dataset_dirs = [p.parent for p in ds_root.glob("*/*/sequences.fasta")]
    if not dataset_dirs:
        raise SystemExit(f"No datasets found under {ds_root}/*/*/sequences.fasta")

    def match_only(spec: str, seg: str) -> bool:
        if not args.only:
            return True
        if "/" not in args.only:
            return False
        a, b = args.only.split("/", 1)
        ok_a = (a == "*" or a == spec)
        ok_b = (b == "*" or b == seg)
        return ok_a and ok_b

    report_rows = []
    for d in sorted(dataset_dirs):
        species = d.parent.name
        segment = d.name
        if not match_only(species, segment):
            continue

        fasta = d / "sequences.fasta"
        meta  = d / "metadata.tsv"
        stats = d / "length_stats.tsv"
        if not (fasta.exists() and meta.exists() and stats.exists()):
            continue

        s = read_length_stats(stats)
        median = int(s["median"])
        p05 = int(s.get("p05", 0))

        min_len = max(p05, int(math.floor(median * args.median_frac)))
        if args.minlen_floor and min_len < args.minlen_floor:
            min_len = args.minlen_floor

        out_dir = filtered_root / species / segment
        out_fasta = out_dir / "sequences.fasta"
        out_meta  = out_dir / "metadata.tsv"
        out_params = out_dir / "filter_params.tsv"

        kept_ids = filter_fasta(fasta, out_fasta, min_len=min_len)
        kept_meta_rows = filter_metadata(meta, out_meta, kept_ids)

        write_filter_params(out_params, {
            "species": species,
            "segment": segment,
            "min_len": min_len,
            "median": median,
            "p05": p05,
            "median_frac": args.median_frac,
            "minlen_floor": args.minlen_floor,
            "kept_fasta": len(kept_ids),
            "kept_meta_rows": kept_meta_rows,
            "source_fasta": str(fasta),
            "source_meta": str(meta),
            "source_stats": str(stats),
        })

        report_rows.append({
            "species": species,
            "segment": segment,
            "min_len": str(min_len),
            "median": str(median),
            "p05": str(p05),
            "kept_fasta": str(len(kept_ids)),
            "kept_meta_rows": str(kept_meta_rows),
            "out_dir": str(out_dir),
        })

    report_path = filtered_root / "filter_report.tsv"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    if report_rows:
        with report_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(report_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(report_rows)

    print(f"Wrote {len(report_rows)} filtered datasets into: {filtered_root}")
    print(f"Report: {report_path}")

if __name__ == "__main__":
    main()