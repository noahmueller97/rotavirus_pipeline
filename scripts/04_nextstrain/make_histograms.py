"""
- Generates sequence length histograms for rotavirus datasets.
- Creates histograms showing sequence length distributions for each species and segment dataset.

Inputs:
- dataset sequences

Outputs:
- histogram plots used for quality assessment
"""

from __future__ import annotations

from pathlib import Path
import argparse
import matplotlib.pyplot as plt


def lengths_from_fasta(path: Path) -> list[int]:
    lengths: list[int] = []
    current_len = 0
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
        if current_len > 0:
            lengths.append(current_len)
    return lengths


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create ONE histogram per dataset folder under results/02_datasets/<species>/<segment>/"
    )
    ap.add_argument(
        "--base-dir",
        type=Path,
        default=Path("results/02_datasets"),
        help="Base datasets directory (default: results/02_datasets)",
    )
    ap.add_argument("--bins", type=int, default=60)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    base: Path = args.base_dir.resolve()
    if not base.exists():
        raise SystemExit(f"Base dir not found: {base}")

    
    fasta_files = sorted(base.glob("*/*/sequences.fasta"))

    if not fasta_files:
        raise SystemExit(f"No datasets found at {base}/<species>/<segment>/sequences.fasta")

    made = 0
    skipped = 0

    for fasta in fasta_files:
        
        try:
            rel = fasta.relative_to(base)
        except ValueError:
            skipped += 1
            continue

        if len(rel.parts) != 3:
            
            skipped += 1
            continue

        species, segment, filename = rel.parts
        if filename != "sequences.fasta":
            skipped += 1
            continue

        out_png = fasta.parent / f"histogram_{species}_{segment}_length_hist.png"
        if out_png.exists() and not args.overwrite:
            skipped += 1
            continue

        lens = lengths_from_fasta(fasta)
        if not lens:
            print(f"[WARN] no sequences in {fasta} -> skipping")
            skipped += 1
            continue

        plt.figure()
        plt.hist(lens, bins=args.bins)
        plt.xlabel("Segment length (nt)")
        plt.ylabel("Count")
        plt.title(f"{species} {segment} length distribution (n={len(lens)})")
        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        plt.close()

        made += 1

    print(f"Done. Created {made} histograms. Skipped {skipped}.")


if __name__ == "__main__":
    main()