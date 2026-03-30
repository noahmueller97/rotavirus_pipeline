"""
- Retrieves coding sequences (CDS) for selected representative accessions.
- Fetches coding sequences for a predefined list of accessions and produces files that can be incorporated into the reference genome configuration.
- Used during reference expansion to identify candidate sequences for improving the reference set.

Inputs:
- accession list embedded in the script or configuration

Outputs:
- top10_reps_cds.tsv
- top10_reps_yaml_snippet.yml

"""
import argparse
from pathlib import Path

from Bio import Entrez, SeqIO

ACCESSIONS = [
    "DQ369973.1",
    "KY634534.1",
    "MH933785.1",
    "PX388204.1",
    "MN029132.1",
    "EU436846.1",
    "OR233779.1",
    "AF516717.1",
    "L18943.1",
    "EU033986.1",
]

def pick_best_cds(record):
    """Return dict for best CDS (by AA length)."""
    best = None
    for feat in record.features:
        if feat.type != "CDS":
            continue
        translation = feat.qualifiers.get("translation", [""])[0]
        if not translation:
            continue

        prot_id = feat.qualifiers.get("protein_id", [""])[0]
        product = feat.qualifiers.get("product", [""])[0]
        gene = feat.qualifiers.get("gene", [""])[0]
        aa_len = len(translation)

        cand = {
            "protein_id": prot_id,
            "product": product,
            "gene": gene,
            "aa_len": aa_len,
        }
        if best is None or cand["aa_len"] > best["aa_len"]:
            best = cand
    return best

def seg_guess_from_product(product: str) -> str:
    """Very light normalization; adjust later if needed."""
    if not product:
        return "unknown"
    p = product.strip().lower()
    
    if "vp1" in p: return "vp1"
    if "vp2" in p: return "vp2"
    if "vp3" in p: return "vp3"
    if "vp4" in p: return "vp4"
    if "vp6" in p: return "vp6"
    if "vp7" in p: return "vp7"
    if "nsp1" in p: return "nsp1"
    if "nsp2" in p: return "nsp2"
    if "nsp3" in p: return "nsp3"
    if "nsp4" in p: return "nsp4"
    if "nsp5" in p or "nsp6" in p: return "nsp5"
    return "unknown"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="NCBI requires a contact email")
    ap.add_argument("--outdir", default="results", help="Output dir (default: results)")
    args = ap.parse_args()

    Entrez.email = args.email

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tsv_path = outdir / "top10_reps_cds.tsv"
    yml_path = outdir / "top10_reps_yaml_snippet.yml"

    rows = []
    for acc in ACCESSIONS:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        best = pick_best_cds(record)
        if best is None:
            rows.append({
                "accession": acc,
                "nucl_len": len(record.seq),
                "protein_id": "",
                "aa_len": 0,
                "product": "",
                "gene": "",
                "segment_guess": "unknown",
            })
            continue

        product = best["product"]
        rows.append({
            "accession": acc,
            "nucl_len": len(record.seq),
            "protein_id": best["protein_id"],
            "aa_len": best["aa_len"],
            "product": product,
            "gene": best["gene"],
            "segment_guess": seg_guess_from_product(product),
        })

    
    with open(tsv_path, "w") as f:
        f.write("accession\tnucl_len\tprotein_id\taa_len\tproduct\tgene\tsegment_guess\n")
        for r in rows:
            f.write(
                f"{r['accession']}\t{r['nucl_len']}\t{r['protein_id']}\t{r['aa_len']}\t"
                f"{r['product']}\t{r['gene']}\t{r['segment_guess']}\n"
            )

   
    with open(yml_path, "w") as f:
        f.write("# Paste these under reference_genomes: (adjust strain/species/virus_type if you want)\n")
        for r in rows:
            if not r["protein_id"]:
                continue
            f.write(
                "  - strain: UNKNOWN_CLUSTER_REP\n"
                "    virus_type: rotavirus\n"
                "    species: unknown\n"
                "    assembly: unknown\n"
                "    segments:\n"
                f"      - name: {r['segment_guess']}\n"
                f"        nucleotide_accession: {r['accession']}\n"
                f"        protein_accession: {r['protein_id']}\n"
                f"        # product: {r['product']}\n"
                "\n"
            )

    print(f"Wrote: {tsv_path}")
    print(f"Wrote: {yml_path}")

if __name__ == "__main__":
    main()