"""
Snakemake pipeline: HLA typing and TCR-HLA merging
"""
import os
from pathlib import Path
from snakemake.exceptions import WorkflowError

# Load configuration file
configfile: "config.yaml"

# Sample list
SAMPLES = config["samples"]

# Path configuration
BAM_DIR = config["bam_dir"]  # /sc/arion/scratch/chos14/NYU/withBAM2025/
HLA_OUTPUT_DIR = config["hla_output_dir"]  # results/hla_typing/
MERGED_OUTPUT_DIR = config["merged_output_dir"]  # results/merged/
SCRIPTS_DIR = "scripts"
THREADS = int(config.get("threads", 16))
CONDA_ENV = "hla-snake"
CONDA_INIT = "/hpc/packages/minerva-rocky9/anaconda3/2025.06/etc/profile.d/conda.sh"
VDJDB_PATH = config["vdjdb_path"]

if not Path(VDJDB_PATH).exists():
    raise WorkflowError(
        f"VDJdb reference not found: {VDJDB_PATH}. "
        "Please download the file and update vdjdb_path in config.yaml."
    )


def get_count_bam(sample_id):
    """
    Count BAM path per sample:
    {bam_dir}/{sample}/outs/per_sample_outs/{sample}/count/sample_alignments.bam
    """
    return os.path.join(
        BAM_DIR,
        sample_id,
        "outs",
        "per_sample_outs",
        sample_id,
        "count",
        "sample_alignments.bam",
    )


def get_tcr_file(sample_id):
    """
    TCR clonotype path per sample:

    {bam_dir}/{sample}/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv
    """
    return os.path.join(
        BAM_DIR,
        sample_id,
        "outs",
        "per_sample_outs",
        sample_id,
        "vdj_t",
        "filtered_contig_annotations.csv",
    )


def get_clonotype_summary(sample_id):
    return os.path.join(MERGED_OUTPUT_DIR, f"{sample_id}_clonotype_summary.tsv")


def get_vdjdb_summary(sample_id):
    return os.path.join(
        MERGED_OUTPUT_DIR, f"{sample_id}_clonotype_summary.vdjdb.tsv"
    )


def get_hla_epitope_summary(sample_id):
    return os.path.join(
        MERGED_OUTPUT_DIR, f"{sample_id}_hla_epitope_summary.tsv"
    )

# All rules
rule all:
    input:
        # Produce the final HLA-epitope summary
        expand(
            os.path.join(MERGED_OUTPUT_DIR, "{sample}_hla_epitope_summary.tsv"),
            sample=SAMPLES,
        )


rule hla_typing:
    """
    Infer HLA genotype from BAM files (via arcasHLA)
    """
    input:
        bam=lambda wildcards: get_count_bam(wildcards.sample)
    output:
        hla_tsv = os.path.join(HLA_OUTPUT_DIR, "{sample}_hla.tsv")
    threads:
        THREADS
    shell:
        """
        # Load cluster modules (optional)
        
        ml proxies || true
        module load micromamba/1.5.3-0 || true
        eval "$(micromamba shell hook --shell bash)" || true
        
        # Activate conda environment
        micromamba activate {CONDA_ENV}
        
        # Reset PYTHONPATH
        unset PYTHONPATH
        
        # Run the arcasHLA helper script
        python {SCRIPTS_DIR}/run_arcas_hla.py \
            --bam {input.bam} \
            --output-dir {HLA_OUTPUT_DIR} \
            --sample-id {wildcards.sample} \
            --threads {threads}
        """


rule merge_tcr_hla:
    """
    Merge TCR clonotype data with HLA typing results
    """
    input:
        tcr_file=lambda wildcards: get_tcr_file(wildcards.sample),
        hla_file=os.path.join(HLA_OUTPUT_DIR, "{sample}_hla.tsv")
    output:
        merged_file = os.path.join(MERGED_OUTPUT_DIR, "{sample}_tcr_hla.tsv")
    shell:
        """
        ml proxies || true
        module load micromamba/1.5.3-0 || true
        eval "$(micromamba shell hook --shell bash)" || true
        micromamba activate {CONDA_ENV}
        python {SCRIPTS_DIR}/merge_tcr_hla.py \
            --tcr-file {input.tcr_file} \
            --hla-file {input.hla_file} \
            --output-file {output.merged_file} \
            --sample-id {wildcards.sample}
        """


rule summarize_clonotypes:
    """
    Aggregate merged TCR-HLA TSV into a clonotype summary
    """
    input:
        merged=os.path.join(MERGED_OUTPUT_DIR, "{sample}_tcr_hla.tsv")
    output:
        summary=os.path.join(MERGED_OUTPUT_DIR, "{sample}_clonotype_summary.tsv")
    shell:
        """
        ml proxies || true
        module load micromamba/1.5.3-0 || true
        eval "$(micromamba shell hook --shell bash)" || true
        micromamba activate {CONDA_ENV}
        python {SCRIPTS_DIR}/summarize_clonotypes.py \
            --input {input.merged} \
            --output {output.summary}
        """


rule annotate_vdjdb:
    """
    Annotate clonotype summary with public epitopes using VDJdb
    """
    input:
        summary=os.path.join(MERGED_OUTPUT_DIR, "{sample}_clonotype_summary.tsv"),
        vdjdb=VDJDB_PATH
    output:
        annotated=os.path.join(
            MERGED_OUTPUT_DIR, "{sample}_clonotype_summary.vdjdb.tsv"
        )
    shell:
        """
        ml proxies || true
        module load micromamba/1.5.3-0 || true
        eval "$(micromamba shell hook --shell bash)" || true
        micromamba activate {CONDA_ENV}
        python {SCRIPTS_DIR}/annotate_vdjdb.py \
            --input {input.summary} \
            --vdjdb {input.vdjdb} \
            --output {output.annotated}
        """


rule summarize_hla_epitopes:
    """
    Summarize at the HLA-epitope level
    """
    input:
        annotated=os.path.join(
            MERGED_OUTPUT_DIR, "{sample}_clonotype_summary.vdjdb.tsv"
        )
    output:
        epitope=os.path.join(MERGED_OUTPUT_DIR, "{sample}_hla_epitope_summary.tsv")
    shell:
        """
        ml proxies || true
        module load micromamba/1.5.3-0 || true
        eval "$(micromamba shell hook --shell bash)" || true
        micromamba activate {CONDA_ENV}
        python {SCRIPTS_DIR}/summarize_hla_epitopes.py \
            --input {input.annotated} \
            --output {output.epitope}
        """

