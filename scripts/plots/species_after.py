"""
- Plots species distribution after dataset filtering.
- Visualizes species composition after filtering and cleaning steps have been applied.

Inputs:
- filtered datasets

Outputs:
- species distribution plot
"""

import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv(
    "results/00_primary/diamond_alignments.tsv",
    sep="\t",
    header=None
)


species = df[1].astype(str).apply(lambda x: x.split("|")[2])

counts = species.value_counts()


plt.figure(figsize=(8,5))

counts.plot(kind="bar")

plt.title("Species distribution after DIAMOND annotation")
plt.xlabel("Rotavirus species")
plt.ylabel("Number of sequences")

plt.tight_layout()
plt.savefig("plots/species_distribution_full.png", dpi=300)
plt.close()



counts_no_A = counts.drop("A", errors="ignore")

plt.figure(figsize=(8,5))

counts_no_A.plot(kind="bar")

plt.title("Species distribution after DIAMOND annotation (excluding A)")
plt.xlabel("Rotavirus species")
plt.ylabel("Number of sequences")

plt.tight_layout()
plt.savefig("plots/species_distribution_noA.png", dpi=300)
plt.close()


print("\nSpecies counts:")
print(counts)