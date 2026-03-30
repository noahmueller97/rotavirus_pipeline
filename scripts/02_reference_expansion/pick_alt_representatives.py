"""
- Selects alternative representative sequences from MMseqs clusters.
- Chooses representative sequences from clusters of unknown sequences to serve as potential new references.
- Used to manually inspect candidate sequences that may improve the reference database.

Inputs:
- cluster TSV from MMseqs
- FASTA file of unknown sequences

Outputs:
- alt_representatives_round2.tsv

"""

import argparse
from collections import defaultdict, Counter
from Bio import SeqIO

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cluster-tsv", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--out", required=True)
    args = p.parse_args()


    members = defaultdict(list)
    with open(args.cluster_tsv) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rep, mem = line.split("\t")[:2]
            members[rep].append(mem)

    
    sizes = {rep: len(mems) for rep, mems in members.items()}
    top_reps = [rep for rep, _ in sorted(sizes.items(), key=lambda x: x[1], reverse=True)[:args.top]]


    wanted_ids = set()
    for rep in top_reps:
        wanted_ids.update(members[rep])

    seq_len = {}
    for rec in SeqIO.parse(args.fasta, "fasta"):
        if rec.id in wanted_ids:
            seq_len[rec.id] = len(rec.seq)


    rows = []
    for rep in top_reps:
        best_id, best_len = None, -1
        for mid in members[rep]:
            L = seq_len.get(mid, -1)
            if L > best_len:
                best_id, best_len = mid, L
        rows.append((rep, sizes[rep], best_id if best_id else "", best_len if best_len >= 0 else ""))


    with open(args.out, "w") as out:
        out.write("rep\tcluster_size\talt_rep_longest\talt_len\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")

    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()