"""
- Plots the number of genome segments per rotavirus strain.
- Shows how many segments are present per strain in the dataset, helping identify incomplete assemblies.

Inputs:
- annotated sequence dataset

Outputs:
- segments-per-strain plot
"""

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "results/00_primary/diamond_alignments.tsv",
    sep="\t",
    header=None
)


strain = df[0]


segment = df[1].astype(str).apply(lambda x: x.split("|")[3])

data = pd.DataFrame({
    "strain": strain,
    "segment": segment
})

segments_per_strain = data.groupby("strain")["segment"].nunique()

plt.figure(figsize=(8,5))

segments_per_strain.value_counts().sort_index().plot(kind="bar")

plt.title("Genome completeness: segments per strain")
plt.xlabel("Number of genome segments")
plt.ylabel("Number of strains")

plt.tight_layout()
plt.savefig("plots/segments_per_strain.png", dpi=300)

plt.show()