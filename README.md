# EpitopeHLA-TCR-snakemake
Single-cell 10x HLA typing and TCR–epitope integration workflow (arcasHLA + VDJdb, Snakemake).

This Snakemake workflow takes **10x Genomics single-cell RNA-seq / V(D)J outputs** and produces a set of tables that link:

* sample-level **HLA genotypes** (via arcasHLA)
* **TCR clonotypes** and per-clonotype aggregates
* **VDJdb** epitope/antigen annotations
* **HLA–epitope summaries** and cohort-level epitope metrics

It is designed for studies of infectious disease and vaccine responses where you want to connect **HLA type ↔ TCR clonotype ↔ known epitope**.

---

## 1. Overview

For each sample in `config.yaml`, the workflow performs:

1. **HLA typing** (arcasHLA)  
   Input: 10x Cell Ranger count BAM  
   Output: `results/hla_typing/{sample_id}_hla.tsv`

2. **TCR–HLA merge**  
   Input: 10x V(D)J contig annotations + HLA typing TSV  
   Output: `results/merged/{sample_id}_tcr_hla.tsv`

3. **Clonotype summarization**  
   Input: merged TCR–HLA TSV  
   Output: `results/merged/{sample_id}_clonotype_summary.tsv`

4. **VDJdb annotation**  
   Input: clonotype summary TSV + VDJdb TSV  
   Output: `results/merged/{sample_id}_clonotype_summary.vdjdb.tsv`

5. **HLA–epitope summary (per sample)**  
   Input: VDJdb-annotated clonotype summary  
   Output: `results/merged/{sample_id}_hla_epitope_summary.tsv`

6. **Cohort-level epitope metrics**  
   Input: all `*_hla_epitope_summary.tsv` + optional `metadata/samples.tsv`  
   Output:  
   * `results/summary/epitope_scores.sample.tsv` – sample-level metrics  
   * `results/summary/epitope_scores.hla.tsv` – sample×HLA-allele-level metrics

These outputs are meant to be consumed by downstream R/Python notebooks for figures
and statistical analysis.

---

## 2. Requirements

- Linux (or any POSIX system where Cell Ranger & arcasHLA can run)
- Conda / Mamba / micromamba
- Snakemake
- arcas-hla
- Python ≥ 3.8
The main conda environment is defined in envs/hla.yaml.

---

## 3. Installation

Clone the repository and create the environment:
```
git clone https://github.com/your-org/your-hla-tcr-pipeline.git
cd your-hla-tcr-pipeline

# Create and activate the environment
conda env create -f envs/hla.yaml
conda activate hla-snake
```
(If you prefer micromamba, replace `conda` with `micromamba`.)

---

## 4. Configuration

### 4.1. `config.yaml`

Edit `config.yaml` to point to your own data. An example is:
```
# List of sample IDs
samples:
  - SAMPLE_001
  - SAMPLE_002

# Root directory with 10x Cell Ranger outputs
bam_dir: "/path/to/cellranger_outputs"

# Output directories
hla_output_dir: "results/hla_typing"
merged_output_dir: "results/merged"

# VDJdb TSV
vdjdb_path: "resources/vdjdb/vdjdb.txt"

# Optional sample metadata table
metadata_tsv: "metadata/samples.tsv"

threads: 16
```

Expected 10x layout for each sample_id:
```
{bam_dir}/{sample_id}/outs/per_sample_outs/{sample_id}/
  ├── count/sample_alignments.bam
  └── vdj_t/filtered_contig_annotations.csv
```

### 4.2. Conda environment

The workflow uses a single environment defined in `envs/hla.yaml`:
```
name: hla-snake
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - arcas-hla
  - snakemake
  - python==3.11
  - pandas
  - biopython
  - pyyaml
```

---

## 5. Running the workflow:

From the project root:
```
conda activate hla-snake

# Dry-run (recommended first)
snakemake -n

# Run locally with conda environments
snakemake --use-conda -j 4
```

To run only up to a specific step, you can request a particular target:
```
# Run up to per-sample HLA typing
snakemake --use-conda -j 4 results/hla_typing/SAMPLE_001_hla.tsv

# Run up to VDJdb-annotated clonotype summary
snakemake --use-conda -j 4 results/merged/SAMPLE_001_clonotype_summary.vdjdb.tsv

# Run the full cohort-level epitope summary
snakemake --use-conda -j 4 results/summary/epitope_scores.sample.tsv
```

On an HPC cluster, you can use a Snakemake profile, for example:
```
snakemake --profile slurm --use-conda -j 100
```

---

## 6. Output
### 6.1. HLA typing results

`results/hla_typing/{sample_id}_hla.tsv`

Tab-delimited, one row per HLA locus:
```
sample_id  locus  allele1      allele2
SAMPLE_001 A      A*01:01      A*02:01
SAMPLE_001 B      B*07:02      B*08:01
SAMPLE_001 C      C*07:01      C*07:02
SAMPLE_001 DRB1   DRB1*03:01   DRB1*04:01
...
```

### 6.2. Merged TCR-HLA table
`results/merged/{sample_id}_tcr_hla.tsv`

Contains all columns from the original 10x V(D)J contig annotations plus per-sample HLA genotype columns

### 6.3. Clonotype summary
