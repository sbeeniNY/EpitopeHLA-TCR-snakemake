#!/usr/bin/env python3
"""
Script to infer HLA genotypes from BAM files via arcasHLA.
"""
import argparse
import subprocess
import sys
import os
import glob
from pathlib import Path


class EmptyFastqError(Exception):
    """Raised when extracted FASTQ files contain no reads."""


def derive_unassigned_bam_path(bam_file):
    """
    Infer the multi/count/unassigned_alignments.bam path from sample_alignments.bam.
    """
    bam_path = Path(bam_file)
    parts = bam_path.parts
    if "outs" not in parts:
        return None
    outs_index = len(parts) - 1 - parts[::-1].index("outs")
    base = Path(*parts[:outs_index + 1])
    fallback_path = base / "multi" / "count" / "unassigned_alignments.bam"
    return str(fallback_path)


def run_arcas_hla(bam_file, output_dir, sample_id, threads=16):
    """
    Run arcasHLA for HLA typing. If no HLA reads are extracted from
    sample_alignments.bam, automatically retry with multi/count/unassigned_alignments.bam.
    """
    os.makedirs(output_dir, exist_ok=True)
    fallback_bam = derive_unassigned_bam_path(bam_file)
    try:
        run_single_bam(
            bam_file=bam_file,
            output_dir=output_dir,
            sample_id=sample_id,
            threads=threads,
            work_suffix=None,
        )
    except EmptyFastqError as err:
        if fallback_bam and os.path.exists(fallback_bam):
            print(
                f"[{sample_id}] No HLA reads detected in {bam_file}. "
                f"Retrying with unassigned BAM: {fallback_bam}"
            )
            run_single_bam(
                bam_file=fallback_bam,
                output_dir=output_dir,
                sample_id=sample_id,
                threads=threads,
                work_suffix="unassigned",
            )
        else:
            raise err


def run_single_bam(bam_file, output_dir, sample_id, threads=16, work_suffix=None):
    """
    Run arcasHLA extract/genotype for a single BAM input.
    """
    work_dir_name = f"{sample_id}_arcas_work"
    if work_suffix:
        work_dir_name = f"{sample_id}_{work_suffix}_arcas_work"
    work_dir = os.path.join(output_dir, work_dir_name)
    os.makedirs(work_dir, exist_ok=True)

    print(f"[{sample_id}] Extracting HLA reads from BAM ({bam_file})...")
    extract_cmd = [
        "arcasHLA", "extract",
        "--single",
        "-t", str(threads),
        "-v",
        "-o", work_dir,
        bam_file
    ]

    result = subprocess.run(extract_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error in arcasHLA extract:\n{result.stderr}")

    bam_path = Path(bam_file)
    bam_stem = bam_path.stem
    bam_name = bam_path.name

    def find_file(patterns):
        for candidate in patterns:
            if os.path.exists(candidate):
                return candidate
        for pattern in patterns:
            matches = glob.glob(pattern, recursive=True)
            if matches:
                return matches[0]
        return None

    fastq_single_candidates = [
        os.path.join(work_dir, f"{bam_stem}.extracted.fq"),
        os.path.join(work_dir, f"{bam_stem}.extracted.fastq"),
        os.path.join(work_dir, f"{sample_id}.extracted.fq"),
        os.path.join(work_dir, f"{sample_id}.extracted.fastq"),
        os.path.join(work_dir, f"{bam_name}.extracted.fq"),
        os.path.join(work_dir, "**", "*.extracted.fq"),
        os.path.join(work_dir, "**", "*.extracted.fastq"),
        os.path.join(work_dir, f"{bam_stem}.extracted.fq.gz"),
        os.path.join(work_dir, f"{bam_stem}.extracted.fastq.gz"),
        os.path.join(work_dir, "**", "*.extracted.fq.gz"),
        os.path.join(work_dir, "**", "*.extracted.fastq.gz"),
    ]

    fastq_single = find_file(fastq_single_candidates)

    extracted_bam_candidates = [
        os.path.join(work_dir, f"{bam_stem}.extracted.bam"),
        os.path.join(work_dir, f"{sample_id}.extracted.bam"),
        os.path.join(work_dir, f"{bam_name}.extracted.bam"),
        os.path.join(work_dir, f"{bam_stem}.bam.extracted.bam"),
        os.path.join(work_dir, "**", "*.extracted.bam"),
    ]
    extracted_bam = find_file(extracted_bam_candidates)

    if fastq_single:
        print(f"[{sample_id}] Genotyping HLA using single FASTQ {os.path.basename(fastq_single)}...")
        genotype_cmd = [
            "arcasHLA", "genotype",
            "--single",
            "--avg", "150",
            "--std", "30",
            "--min_count", "20",
            "-g", "A,B,C,DQA1,DQB1,DRB1",
            "-t", str(threads),
            "-v",
            "-o", work_dir,
            fastq_single,
        ]
    elif extracted_bam:
        print(f"[{sample_id}] Genotyping HLA using {extracted_bam}...")
        genotype_cmd = [
            "arcasHLA", "genotype",
            "-t", str(threads),
            "-v",
            "-o", work_dir,
            extracted_bam
        ]
    else:
        raise RuntimeError(
            f"No extracted FASTQ or BAM files found in {work_dir}. "
            "Checked for patterns like *.extracted.1.fq, *.extracted.2.fq, *.extracted.bam"
        )

    result = subprocess.run(genotype_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        stderr = result.stderr or ""
        if "FASTQ files are empty" in stderr:
            raise EmptyFastqError(stderr)
        raise RuntimeError(f"Error in arcasHLA genotype:\n{stderr}")

    output_file = os.path.join(output_dir, f"{sample_id}_hla.tsv")
    possible_genotype_files = [
        os.path.join(work_dir, f"{bam_stem}.genotype.json"),
        os.path.join(work_dir, f"{bam_name}.genotype.json"),
        os.path.join(work_dir, f"{sample_id}.genotype.json"),
    ]
    if extracted_bam:
        possible_genotype_files.append(
            os.path.join(work_dir, f"{Path(extracted_bam).stem}.genotype.json")
        )

    genotype_file = None
    for f in possible_genotype_files:
        if os.path.exists(f):
            genotype_file = f
            break

    if genotype_file is None:
        found_files = glob.glob(os.path.join(work_dir, "**", "*.genotype.json"), recursive=True)
        if found_files:
            genotype_file = found_files[0]

    if genotype_file and os.path.exists(genotype_file):
        convert_json_to_tsv(genotype_file, output_file, sample_id)
        print(f"[{sample_id}] HLA typing completed: {output_file}")
    else:
        raise RuntimeError(
            f"Genotype file not found in {work_dir}. Looked for files matching *.genotype.json"
        )


def convert_json_to_tsv(json_file, tsv_file, sample_id):
    """
    Convert arcasHLA JSON outputs into TSV format at the 2-field level (e.g., A*01:01).
    """
    import json
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # List of HLA loci
    loci = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    
    results = []
    for locus in loci:
        if locus in data:
            alleles = data[locus]
            # Convert to a 2-field level (e.g., A*01:01:01:01 -> A*01:01)
            if isinstance(alleles, list) and len(alleles) >= 2:
                allele1 = alleles[0].split(':')[:2]  # Keep only the first two fields
                allele2 = alleles[1].split(':')[:2]
                allele1_str = f"{locus}*{':'.join(allele1)}"
                allele2_str = f"{locus}*{':'.join(allele2)}"
                results.append({
                    'sample_id': sample_id,
                    'locus': locus,
                    'allele1': allele1_str,
                    'allele2': allele2_str
                })
    
    # Save as TSV
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(tsv_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run arcasHLA for HLA typing")
    parser.add_argument("--bam", required=True, help="Input BAM file path")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    
    args = parser.parse_args()
    
    try:
        run_arcas_hla(
            bam_file=args.bam,
            output_dir=args.output_dir,
            sample_id=args.sample_id,
            threads=args.threads
        )
    except EmptyFastqError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

