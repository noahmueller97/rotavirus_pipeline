"""
- Annotates rotavirus sequences using DIAMOND alignment results.
- Maps DIAMOND protein hits to rotavirus genome segments and adds segment annotations to each sequence record.

Inputs:
- genbank.ndjson
- diamond_alignments.tsv

Outputs:
- annotated_genbank.ndjson
"""

import pandas as pd

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate NCBI ndjson with rotavirus segment/species info from diamond tsv")
    parser.add_argument("--ndjson", type=str, required=True, help="Input ndjson file with GenBank records")
    parser.add_argument("--diamond-tsv", type=str, required=True, help="Diamond tsv file with annotations")
    parser.add_argument("--pident-threshold", type=float, default=70, help="Minimum percent identity to consider annotation valid")
    parser.add_argument("--output", dest="output_ndjson", type=str, required=True, help="Output ndjson file")
    args = parser.parse_args()

   
    df = pd.read_csv(
        args.diamond_tsv,
        header=None,
        names=['accession', 'match', 'pident', 'length', 'mismatch',
               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'],
        sep="\t"
    )
    df.set_index('accession', inplace=True)

    segment_lookup = {}
    for _, row in df.iterrows():
        accession = row.name

        
        
        parts = row['match'].split('|')
        if len(parts) < 4:
            continue

        ref_nuc_acc = parts[0]
        virus_type = parts[1] if len(parts) > 1 else None
        species = parts[2] if len(parts) > 2 else None
        segment = parts[3] if len(parts) > 3 else None

        segment_lookup[accession] = {
            'new_segment': segment,
            'new_species': species,
            'match': ref_nuc_acc,
            'virus_type': virus_type,
            'pident': row['pident'],
        }

    
    from augur.io.json import load_ndjson
    import json
    from collections import defaultdict

    missing_records = []
    mapping_counts = defaultdict(int)

    with open(args.output_ndjson, 'w') as out_f, open(args.ndjson, 'r') as f:
        for record in load_ndjson(f):
            seq = record.get('sequence', '')
            sequence_length = len(seq)
            n_count = seq.upper().count('N')
            good_sequence = 1 - (n_count / sequence_length if sequence_length > 0 else 0)

            accession = record.get('Accession')
            if accession in segment_lookup and (segment_lookup[accession]['pident'] / max(good_sequence, 1e-6)) > args.pident_threshold:
                record.update(segment_lookup[accession])
            else:
                missing_records.append(accession)

            mapping_counts[(record.get('Segment', 'unknown'), record.get('new_segment', 'unknown'))] += 1
            out_f.write(json.dumps(record) + "\n")

    print(f"Missing records for {len(missing_records)} accessions")
    for (old_seg, new_seg), count in mapping_counts.items():
        print(f"Segment mapping {old_seg} -> {new_seg}: {count} records")