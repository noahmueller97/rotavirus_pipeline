"""
- Exports a rotavirus dataset in standardized format.
- Processes intermediate data structures and exports them into
a standardized dataset format used by downstream pipeline
steps.

Inputs:
- intermediate dataset structures

Outputs:
- dataset files used by the Nextstrain preparation steps
"""

import json
import argparse
from pathlib import Path

def norm(x):
    if x is None:
        return ""
    return str(x).strip()

def is_rva(r):
    
    ns = norm(r.get("new_species")).upper()
    if ns:
        return ns == "A"
    
    return norm(r.get("Isolate Lineage")).startswith("RVA/")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--ndjson", required=True)
    p.add_argument("--species", default="A")     
    p.add_argument("--segment", required=True)    
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    seg = args.segment.lower()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_path = outdir / "sequences.fasta"
    meta_path  = outdir / "metadata.tsv"

    n_total = n_match = n_written = 0

    with open(fasta_path, "w") as fasta, open(meta_path, "w") as meta, open(args.ndjson) as f:
        meta.write("strain\tdate\tcountry\tregion\thost\tsegment\tlineage\tnew_species\tpident\n")

        for line in f:
            n_total += 1
            r = json.loads(line)

            if args.species.upper() == "A" and not is_rva(r):
                continue

            if norm(r.get("new_segment")).lower() != seg:
                continue

            n_match += 1

            acc = norm(r.get("Accession"))
            seq = norm(r.get("sequence")).replace("\n", "")
            if not acc or not seq:
                continue

            date = norm(r.get("Isolate Collection date"))
            loc  = norm(r.get("Geographic Location"))
            reg  = norm(r.get("Geographic Region"))
            host = norm(r.get("Host Name"))
            lineage = norm(r.get("Isolate Lineage"))
            new_species = norm(r.get("new_species"))
            pident = norm(r.get("pident"))

            country = loc.split(":")[0].strip() if loc else ""

            fasta.write(f">{acc}\n{seq}\n")
            meta.write(f"{acc}\t{date}\t{country}\t{reg}\t{host}\t{seg}\t{lineage}\t{new_species}\t{pident}\n")
            n_written += 1

    print("Total records:", n_total)
    print("Matched (species+segment):", n_match)
    print("Written (has accession+sequence):", n_written)
    print("FASTA:", fasta_path)
    print("META :", meta_path)