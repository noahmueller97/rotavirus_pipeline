"""
- Exports sequences for a single genome segment.
- Extracts sequences belonging to a specific rotavirus genome segment and writes them to FASTA and metadata files.

Inputs:
- annotated sequence dataset

Outputs:
- segment-specific FASTA
- segment-specific metadata
"""

import json
from pathlib import Path

NDJSON  = Path("results/annotated_genbank.ndjson")
SEGMENT = "nsp5"
OUTDIR  = Path(f"results/trees/{SEGMENT}")

def norm(x):
    return str(x).strip().lower() if x else ""

OUTDIR.mkdir(parents=True, exist_ok=True)
fasta_out = OUTDIR / f"{SEGMENT}.fasta"

n = 0
with NDJSON.open() as f, fasta_out.open("w") as out:
    for line in f:
        r = json.loads(line)
        seg = norm(r.get("new_segment"))
        if seg != SEGMENT:
            continue
        acc = (r.get("Accession") or "").strip()
        seq = (r.get("sequence") or "").replace("\n","").strip()
        if not acc or not seq:
            continue
        out.write(f">{acc}\n{seq}\n")
        n += 1

print("wrote", fasta_out, "n=", n)