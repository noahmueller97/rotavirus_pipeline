# Rotavirus Phylogenetic Analysis Pipeline

This repository contains a reproducible Snakemake pipeline for downloading, annotating, filtering, and phylogenetically analysing rotavirus genome sequences. The workflow retrieves rotavirus sequences from NCBI, annotates gene segments using DIAMOND, filters and subsamples datasets, and reconstructs phylogenetic trees using the Nextstrain toolkit.

The pipeline was developed to analyse the evolutionary relationships of rotavirus strains and to generate interactive phylogenetic visualizations.


Repository Structure

config/      configuration files used by the pipeline
scripts/     Python scripts used in different workflow stages
workflow/    Snakemake workflow definitions


Scripts

The scripts directory contains helper scripts grouped by pipeline stage.

scripts/
 01_annotation
 02_reference_expansion
 03_dataset_building
 04_nextstrain
 plots


01_annotation
Scripts for annotating sequences and calculating basic statistics.

Examples include:
- annotate_sequences.py
- basic_stats.py
- fetch_reference_genomes.py


02_reference_expansion
Identifies sequences that could not be classified during DIAMOND alignment. These sequences are clustered using MMseqs2 and representative sequences are selected to expand the reference database.

This step was primarily used during development to improve reference coverage.


03_dataset_building
Builds cleaned and filtered datasets for downstream phylogenetic analysis.


04_nextstrain
Scripts used in the Nextstrain analysis stage, including filtering, statistics, and subsampling of datasets.


plots
Scripts used to generate summary plots for exploratory analysis and figures.


Workflow Overview

The pipeline is organized into multiple Snakemake workflows located in the workflow directory.

workflow/
 Snakefile_01_data
 Snakefile_02_reference
 Snakefile_03_filter
 Snakefile_04_cleaning
 Snakefile_05_nextstrain

Each Snakefile represents a stage in the analysis pipeline.


Pipeline Stages


1. Data acquisition and annotation

Snakefile_01_data

This step:
- downloads rotavirus sequences from NCBI
- extracts FASTA sequences and metadata
- builds a DIAMOND reference database
- aligns sequences against reference proteins
- annotates genome segments
- generates basic statistics

Output examples:

results/00_primary/
 diamond_alignments.tsv
 segment_comparison.tsv
 annotated_genbank.ndjson

Run:
snakemake -s workflow/Snakefile_01_data -j 4


2. Reference expansion (optional)

Snakefile_02_reference

This stage identifies sequences that could not be classified during DIAMOND annotation.

Steps include:
1. extracting unknown sequences
2. clustering them with MMseqs2
3. identifying representative sequences
4. retrieving candidate coding sequences

The resulting sequences can be used to extend the reference database.

Run:
snakemake -s workflow/Snakefile_02_reference -j 4

Outputs are written to:
results/01_unknowns/

This step is optional once the reference set has been finalized.


3. Dataset filtering

Snakefile_03_filter

This stage prepares datasets for phylogenetic analysis by:
- filtering sequences based on length
- ensuring sequences contain the expected gene segments
- generating clean datasets for each species and segment

Outputs:

results/03_nextstrain/01_filtered/
results/03_nextstrain/02_subsampled/

Run:
snakemake -s workflow/Snakefile_03_filter -j 4


4. Dataset cleaning

Snakefile_04_cleaning

Performs additional quality control and prepares metadata required for downstream phylogenetic analysis.

Run:
snakemake -s workflow/Snakefile_04_cleaning -j 4


5. Phylogenetic reconstruction and visualization

Snakefile_05_nextstrain

This stage performs the phylogenetic analysis using the Nextstrain toolkit.

Steps include:
1. sequence alignment using augur align
2. phylogenetic tree inference using augur tree
3. tree refinement using augur refine
4. exporting interactive visualization datasets

Outputs:

results/03_nextstrain/03_aligned/
results/03_nextstrain/04_tree/
results/03_nextstrain/05_auspice/

The final output is a Nextstrain-compatible JSON dataset.

Example:
results/03_nextstrain/05_auspice/A/vp1/auspice.json

Run:
snakemake -s workflow/Snakefile_05_nextstrain -j 4

To generate a single example tree:
snakemake -s workflow/Snakefile_05_nextstrain results/03_nextstrain/05_auspice/A/vp1/auspice.json -j 4


Visualizing the phylogenetic trees

The resulting datasets can be visualized using the Nextstrain Auspice interface.

Run:
npx auspice view --datasetDir results/03_nextstrain/05_auspice/<Strain>/<Segment>

Then open:
http://localhost:4000

The interface allows interactive exploration of the phylogenetic trees, including filtering by geographic region, host species, and other metadata attributes.


Dependencies

The pipeline relies on the following tools:

- Snakemake
- Python
- DIAMOND
- MMseqs2
- Nextstrain / Augur
- Biopython


Notes

- The pipeline is modular and can be executed stage-by-stage.
- Intermediate outputs are stored in the results directory.
- The reference expansion stage is optional once a stable reference database has been established.