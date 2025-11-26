# HLA Typing Pipeline Guide

This workflow estimates HLA genotypes from 10x single-cell RNA-seq BAM files and merges those calls with TCR clonotype data.

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Configuration](#configuration)
5. [Running](#running)
6. [Outputs](#outputs)
7. [Troubleshooting](#troubleshooting)
8. [Standalone Scripts](#standalone-scripts)
9. [References](#references)

## Overview

The Snakemake pipeline runs the following steps:

1. **HLA typing** – Use arcasHLA on each 10x Cell Ranger count BAM (`{bam_dir}/{sample_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam`).  
   Output: `results/hla_typing/{sample_id}_hla.tsv`

2. **TCR–HLA merge** – Combine the TCR clonotype table (`{bam_dir}/{sample_id}/outs/per_sample_outs/{sample_id}/vdj_t/filtered_contig_annotations.csv`) with the HLA typing results.  
   Output: `results/merged/{sample_id}_tcr_hla.tsv`

3. **Clonotype summary** – Aggregate `sample_id + chain + clonotype_id` rows via `scripts/summarize_clonotypes.py` (counts cells/contigs, picks representative CDR3/V/J, sums reads/UMIs, carries HLA QC columns).  
   Output: `results/merged/{sample_id}_clonotype_summary.tsv`

4. **VDJdb annotation** – Match public epitope/antigen information with `scripts/annotate_vdjdb.py`.  
   Input: clonotype summary + VDJdb TSV (`resources/vdjdb/vdjdb.txt` or the `vdjdb_path` from `config.yaml`).  
   Output: `results/merged/{sample_id}_clonotype_summary.vdjdb.tsv`

5. **HLA–epitope summary** – Collapse annotated clonotypes by (sample, HLA allele, antigen) with `scripts/summarize_hla_epitopes.py`.  
   Output: `results/merged/{sample_id}_hla_epitope_summary.tsv`

## Requirements

- Linux or another POSIX environment
- Conda / Mamba
- Snakemake
- arcasHLA (installed through Conda)
- Python 3.8+

## Installation

### 1. Create the Conda environment

```bash
cd HLA_typing
conda env create -f envs/hla.yaml
conda activate hla-snake
# micromamba users:
# micromamba create -f envs/hla.yaml
```

### 2. Prepare arcasHLA reference data (one-time)

```bash
conda activate arcas-hla
unset PYTHONPATH
arcasHLA reference
```

If the reference bundle is missing, download the latest release from [IMGTHLA](https://github.com/ANHIG/IMGTHLA) and install it manually for arcasHLA.

## Configuration

### 1. Edit `config.yaml`

Adjust the configuration to point at your resources:

```yaml
# Sample list
samples:
  - P029_D3
  - P029_D7_HIMC

# Root directory containing 10x outputs
bam_dir: "/sc/arion/scratch/chos14/NYU/withBAM2025"  # {sample}/outs/per_sample_outs/... root

# Output directories
hla_output_dir: "/sc/arion/projects/scMultiscale/chos14/NYU/05_HLA/"
merged_output_dir: "/sc/arion/projects/scMultiscale/chos14/NYU/05_HLA/merged"

# Optional VDJdb TSV
vdjdb_path: "resources/vdjdb/vdjdb.txt"

# Runtime options
threads: 16
```

### 2. Verify directory layout

Ensure the repository and input files follow this structure:

```
HLA_typing/
├── Snakefile
├── config.yaml
├── scripts/
│   ├── run_arcas_hla.py
│   └── merge_tcr_hla.py
├── envs/
│   └── hla.yaml
└── PIPELINE.md

Input root: /sc/arion/scratch/chos14/NYU/withBAM2025/
  ├── P029_D3/
  │   └── outs/per_sample_outs/P029_D3/
  │       ├── count/sample_alignments.bam
  │       └── vdj_t/filtered_contig_annotations.csv
  └── P029_D7_HIMC/outs/per_sample_outs/P029_D7_HIMC/{count,vdj_t}/...
```

## Running

### Full pipeline

```bash
cd HLA_typing
conda activate hla-snake

# Dry-run
snakemake -n

# Execute
snakemake -j 4
```

### Specific samples

```bash
snakemake results/merged/sample1_tcr_hla.tsv -j 1
```

### Let Snakemake manage Conda environments

```bash
snakemake --use-conda -j 1
```

### Cluster execution

On SLURM or similar schedulers, use a profile:

```bash
snakemake --profile slurm -j 100
```

## Outputs

### 1. HLA typing results

`results/hla_typing/{sample_id}_hla.tsv`

Format:
```
sample_id  locus  allele1      allele2
sample1    A      A*01:01      A*02:01
sample1    B      B*07:02      B*08:01
sample1    C      C*07:01      C*07:02
sample1    DRB1   DRB1*03:01   DRB1*04:01
...
```

- Reported at the **2-field level** (e.g., `A*01:01`)
- Main loci: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1

### 2. Merged TCR–HLA table

`results/merged/{sample_id}_tcr_hla.tsv`

- Contains every column from the original TCR clonotype table
- Adds HLA genotype columns:
  - `sample_id`
  - `HLA_A_allele1`, `HLA_A_allele2`
  - `HLA_B_allele1`, `HLA_B_allele2`
  - `HLA_C_allele1`, `HLA_C_allele2`
  - `HLA_DRB1_allele1`, `HLA_DRB1_allele2`
  - ...and equivalent columns for other loci

### 3. Clonotype summary

`results/merged/{sample_id}_clonotype_summary.tsv`

- Sample columns: `sample_id`, `clonotype_id`, `chain`, `n_cells`, `n_contigs`, `cdr3_aa`, `v_gene`, `total_reads`, `total_umis`, `HLA_*_allele{1,2}`
- Falls back to `raw_clonotype_id` when `exact_subclonotype_id` is missing

### 4. VDJdb annotation results

`results/merged/{sample_id}_clonotype_summary.vdjdb.tsv`

- Adds `vdjdb_match`, `vdjdb_n_matches`, `vdjdb_antigen`, `vdjdb_antigen_species`, `vdjdb_mhc`, `vdjdb_reference`
- Adjust column mappings with the CLI options in `annotate_vdjdb.py` if your VDJdb file uses different headers

### 5. HLA–epitope summary

`results/merged/{sample_id}_hla_epitope_summary.tsv`

- Provides `n_clonotypes` and `n_cells` totals per (sample, HLA allele, VDJdb antigen) combination
- Helpful for quickly assessing public clonotype expansion across HLA backgrounds

## Troubleshooting

### 1. Missing arcasHLA reference

Error message:
```
FileNotFoundError: ... hla.dat
```

Resolution:
- Download the latest [IMGTHLA](https://github.com/ANHIG/IMGTHLA) release
- Place the contents inside the arcasHLA installation directory

### 2. BAM file not found

- Confirm the `bam_dir` value inside `config.yaml`
- Verify the BAM filename matches `{sample_id}.bam`

### 3. TCR file not found

- Check for `{bam_dir}/{sample_id}/outs/per_sample_outs/{sample_id}/vdj_t/filtered_contig_annotations.csv`
- Ensure the Cell Ranger V(D)J pipeline completed successfully

### 4. Conda activation fails

If your server requires environment modules, update the `module load` commands inside the `hla_typing` rule in the `Snakefile`.

### 5. Permission errors

Confirm write access to the output directories:

```bash
mkdir -p results/hla_typing results/merged
```

## Standalone Scripts

### Run only HLA typing

```bash
python scripts/run_arcas_hla.py \
    --bam /sc/arion/scratch/chos14/NYU/withBAM2025/sample1/outs/per_sample_outs/sample1/count/sample_alignments.bam \
    --output-dir /sc/arion/projects/scMultiscale/chos14/NYU/05_HLA/ \
    --sample-id sample1 \
    --threads 16
```

### Run only the TCR–HLA merge

```bash
python scripts/merge_tcr_hla.py \
    --tcr-file /sc/arion/scratch/chos14/NYU/withBAM2025/sample1/outs/per_sample_outs/sample1/vdj_t/filtered_contig_annotations.csv \
    --hla-file /sc/arion/projects/scMultiscale/chos14/NYU/05_HLA/sample1_hla.tsv \
    --output-file /sc/arion/projects/scMultiscale/chos14/NYU/05_HLA/merged/sample1_tcr_hla.tsv \
    --sample-id sample1
```

## References

- [arcasHLA GitHub](https://github.com/RabadanLab/arcasHLA)
- [Snakemake documentation](https://snakemake.readthedocs.io/)
- [IMGTHLA](https://github.com/ANHIG/IMGTHLA)
