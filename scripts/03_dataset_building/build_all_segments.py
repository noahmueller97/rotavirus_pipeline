"""
- Builds per-species and per-segment rotavirus datasets.
- Constructs dataset directories for each species and genome segment by splitting the annotated NDJSON dataset.
- This script prepares the structured datasets used for subsequent filtering and phylogenetic analysis.

Inputs:
- annotated_genbank.ndjson
- diamond_alignments.tsv

Outputs:
- datasets organized by species and segment
- results/02_datasets/

"""

import json
from pathlib import Path
from collections import defaultdict
import re

NDJSON  = Path("results/00_primary/annotated_genbank.ndjson")
DIAMOND = Path("results/00_primary/diamond_alignments.tsv")
OUTDIR  = Path("results/02_datasets")

def norm(x):
    if x is None:
        return ""
    return str(x).strip()

def norm_lower(x):
    return norm(x).lower()

def parse_species_from_lineage(lineage: str) -> str:
    """
    Typical: RVA/Human-wt/CHN/... -> species A
             RVB/... -> B, etc.
    Returns e.g. "A", "B", "C" or "" if unknown.
    """
    if not lineage:
        return ""
    m = re.match(r"^\s*RV([A-Z])\b", lineage.strip())
    return m.group(1) if m else ""


coords = {}
with DIAMOND.open() as f:
    for line in f:
        if not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8:
            continue
        acc = cols[0]
        qstart = int(cols[6]) - 1  
        qend   = int(cols[7])      
        coords[acc] = (qstart, qend)

OUTDIR.mkdir(parents=True, exist_ok=True)

writers = {}               
seen = defaultdict(set)    

total = 0
written = 0
skipped_no_coords = 0
skipped_no_segment = 0
skipped_no_species = 0
skipped_short = 0

with NDJSON.open() as f:
    for line in f:
        if not line.strip():
            continue
        total += 1
        r = json.loads(line)

        acc = norm(r.get("Accession"))
        if not acc:
            continue

        seg = norm_lower(r.get("new_segment"))
        if not seg or seg == "unknown":
            skipped_no_segment += 1
            continue

        
        species = norm(r.get("new_species"))
        if not species or species.lower() in ("unknown", "null", "none"):
            species = parse_species_from_lineage(norm(r.get("Isolate Lineage")))

        species = species.upper().strip()
        if not species:
            skipped_no_species += 1
            continue

        full_seq = norm(r.get("sequence")).replace("\n", "")
        if not full_seq:
            continue

        if acc not in coords:
            skipped_no_coords += 1
            continue

        start, end = coords[acc]
        if start < 0 or end > len(full_seq) or start >= end:
            skipped_no_coords += 1
            continue

        seq = full_seq[start:end]
        if len(seq) < 200:
            skipped_short += 1
            continue

        key = (species, seg)
        if acc in seen[key]:
            continue
        seen[key].add(acc)

       
        segdir = OUTDIR / species / seg
        segdir.mkdir(parents=True, exist_ok=True)

        fasta_path = segdir / "sequences.fasta"
        meta_path  = segdir / "metadata.tsv"

        if key not in writers:
            fasta = open(fasta_path, "w")
            meta  = open(meta_path, "w")
            meta.write("strain\tdate\tcountry\tregion\thost\tspecies\tsegment\n")
            writers[key] = (fasta, meta)

        fasta, meta = writers[key]

        date   = norm(r.get("Isolate Collection date"))

        loc    = norm(r.get("Geographic Location"))
        region = norm(r.get("Geographic Region"))
        host   = norm(r.get("Host Name"))

        
        country = loc.split(":")[0].strip() if loc else ""
        if not country:
            country = norm(r.get("Submitter Country"))

        fasta.write(f">{acc}\n{seq}\n")
        meta.write(f"{acc}\t{date}\t{country}\t{region}\t{host}\t{species}\t{seg}\n")
        written += 1

for fasta, meta in writers.values():
    fasta.close()
    meta.close()

print("Done.")
print("Total records read:", total)
print("Written:", written)
print("Skipped (no segment / unknown):", skipped_no_segment)
print("Skipped (no species):", skipped_no_species)
print("Skipped (no diamond coords / bad slice):", skipped_no_coords)
print("Skipped (too short):", skipped_short)
print("Output root:", OUTDIR)