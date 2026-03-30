"""
- Generates host distribution plots for rotavirus sequences.
- Visualizes the distribution of host species across the dataset.

Inputs:
- annotated sequence dataset

Outputs:
- host distribution plot
"""

import pandas as pd
import matplotlib.pyplot as plt
import glob

files = glob.glob("results/03_nextstrain/01_filtered/*/*/metadata.tsv")

dfs = [pd.read_csv(f, sep="\t") for f in files]

data = pd.concat(dfs)

hosts = data["host"].value_counts()

plt.figure(figsize=(8,5))

hosts.head(10).plot(kind="bar")

plt.title("Host distribution")
plt.xlabel("Host species")
plt.ylabel("Number of sequences")

plt.tight_layout()
plt.savefig("plots/host_distribution.png", dpi=300)

plt.show()