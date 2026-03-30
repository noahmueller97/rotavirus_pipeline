"""
- Fetches rotavirus reference proteins from configured accessions.
- Downloads reference genome/protein sequences defined in config/reference_genomes.yml and writes them into a FASTA file used to build the DIAMOND database.
- Reference proteins are later used to annotate genome segments via DIAMOND alignment.

Inputs:
- reference_genomes.yml (configuration file)

Outputs:
- reference_proteins.fasta

"""

import yaml
from Bio import Entrez, SeqIO, Seq

Entrez.email = "noah.mueller@stud.unibas.ch"

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fetch reference genomes from NCBI")
    parser.add_argument("--references", type=str, required=True, help="yaml config file with accessions")
    parser.add_argument("--out-genomes", type=str, default="data/reference_genomes.fasta")
    parser.add_argument("--out-proteins", type=str, default="data/reference_proteins.fasta")
    args = parser.parse_args()

    with open(args.references, 'r') as f:
        references = yaml.safe_load(f)

    seqs = []
    aa_seqs = []
    for ref in references['reference_genomes']:
        strain = ref['strain']
        virus_type = ref.get('virus_type', 'rotavirus')
        species = ref.get('species', None)

        for segment in ref['segments']:
            accession = segment.get('nucleotide_accession') or segment.get('accession')
            segment_name = segment['name']

            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            
            id_parts = [accession, virus_type, species, segment_name]
            record.id = "|".join([p for p in id_parts if p])
            record.description = f"strain={strain}"
            seqs.append(record)

            for feature in record.features:
                if feature.type == "CDS":
                    aa_seq = feature.qualifiers.get("translation", [""])[0]
                    if aa_seq:
                        cds_name = feature.qualifiers.get('protein_id', ['unknown'])[0]
                        aa_id_parts = [accession, virus_type, species, segment_name, cds_name]
                        aa_record = SeqIO.SeqRecord(
                            seq=Seq.Seq(aa_seq),
                            id="|".join([p for p in aa_id_parts if p]),
                            description=f"strain={strain}"
                        )
                        aa_seqs.append(aa_record)

    SeqIO.write(seqs, args.out_genomes, "fasta")
    SeqIO.write(aa_seqs, args.out_proteins, "fasta")