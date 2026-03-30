"""
- Exports unknown rotavirus sequences for downstream analysis.
- Extracts sequences that could not be assigned to known genome segments during the annotation step.
- Unknown sequences are clustered to identify potential missing reference segments.

Inputs:
- annotated_genbank.ndjson

Outputs:
- FASTA file containing unknown sequences

"""

import json

INFILE = "/Users/noahmueller/virus-nextstrain/rotavirus/results/annotated_genbank.ndjson"
OUTFILE = "/Users/noahmueller/virus-nextstrain/rotavirus/results/diamond_unknowns.tsv"


all_keys = set()
records = []

with open(INFILE) as f:
    for line in f:
        r = json.loads(line)

        
        if r.get("new_segment"):
            continue

        records.append(r)
        all_keys.update(r.keys())


cols = sorted(all_keys)

with open(OUTFILE, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in records:
        out.write(
            "\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols)
            + "\n"
        )

print(f"Wrote {len(records)} DIAMOND-unknown records to {OUTFILE}")