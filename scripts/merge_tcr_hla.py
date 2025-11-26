#!/usr/bin/env python3
"""
Script for merging TCR clonotype tables with HLA typing results.
"""
import argparse
import pandas as pd
import sys
import os


def load_hla_results(hla_file):
    """
    Load HLA typing results and aggregate them per donor.

    Returns
    -------
    dict
        Dictionary shaped as {locus: [allele1, allele2]}.
    """
    hla_df = pd.read_csv(hla_file, sep='\t')
    
    # Group by sample_id (a donor can have multiple samples)
    hla_dict = {}
    for _, row in hla_df.iterrows():
        locus = row['locus']
        if locus not in hla_dict:
            hla_dict[locus] = []
        hla_dict[locus].extend([row['allele1'], row['allele2']])
    
    # Deduplicate and sort alleles for each locus
    for locus in hla_dict:
        hla_dict[locus] = sorted(list(set(hla_dict[locus])))
    
    return hla_dict


def load_tcr_data(tcr_file):
    """
    Load a 10x Genomics TCR clonotype table.
    """
    tcr_df = pd.read_csv(tcr_file)
    return tcr_df


def merge_tcr_hla(tcr_file, hla_file, output_file, sample_id):
    """
    Merge TCR data with HLA results.

    Parameters
    ----------
    tcr_file : str
        Path to the TCR clonotype file ({bam_dir}/{sample_id}/outs/per_sample_outs/{sample_id}/vdj_t/filtered_contig_annotations.csv)
    hla_file : str
        Path to the HLA typing result (results/hla_typing/{sample_id}_hla.tsv)
    output_file : str
        Destination for the merged file (results/merged/{sample_id}_tcr_hla.tsv)
    sample_id : str
        Sample identifier
    """
    # Verify the input files exist
    if not os.path.exists(tcr_file):
        print(f"Error: TCR file not found: {tcr_file}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(hla_file):
        print(f"Error: HLA file not found: {hla_file}", file=sys.stderr)
        sys.exit(1)
    
    # Load data
    print(f"[{sample_id}] Loading TCR data...")
    tcr_df = load_tcr_data(tcr_file)
    
    print(f"[{sample_id}] Loading HLA data...")
    hla_dict = load_hla_results(hla_file)
    
    # Attach sample_id to the TCR dataframe
    tcr_df['sample_id'] = sample_id
    
    # Append HLA information as columns
    hla_loci = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    for locus in hla_loci:
        if locus in hla_dict:
            alleles = hla_dict[locus]
            if len(alleles) >= 2:
                tcr_df[f'HLA_{locus}_allele1'] = alleles[0]
                tcr_df[f'HLA_{locus}_allele2'] = alleles[1]
            elif len(alleles) == 1:
                tcr_df[f'HLA_{locus}_allele1'] = alleles[0]
                tcr_df[f'HLA_{locus}_allele2'] = alleles[0]
        else:
            tcr_df[f'HLA_{locus}_allele1'] = None
            tcr_df[f'HLA_{locus}_allele2'] = None
    
    # Create the output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Persist merged results
    tcr_df.to_csv(output_file, sep='\t', index=False)
    print(f"[{sample_id}] Merged TCR-HLA data saved: {output_file}")
    print(f"  - TCR records: {len(tcr_df)}")
    print(f"  - HLA loci: {len([l for l in hla_loci if l in hla_dict])}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge TCR clonotype and HLA typing results")
    parser.add_argument("--tcr-file", required=True, help="TCR clonotype file path")
    parser.add_argument("--hla-file", required=True, help="HLA typing result file path")
    parser.add_argument("--output-file", required=True, help="Output merged file path")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    
    args = parser.parse_args()
    
    merge_tcr_hla(
        tcr_file=args.tcr_file,
        hla_file=args.hla_file,
        output_file=args.output_file,
        sample_id=args.sample_id
    )

