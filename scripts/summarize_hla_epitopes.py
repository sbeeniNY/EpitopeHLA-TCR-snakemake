#!/usr/bin/env python3
"""
Summarize VDJdb-annotated clonotypes by HLA allele and antigen.

Input:
    TSV produced by annotate_vdjdb.py. Requires the columns:
        sample_id, n_cells, vdjdb_match, vdjdb_antigen,
        vdjdb_antigen_species, vdjdb_entries (JSON list),
        HLA_*_allele1/2.

Output:
    A TSV with columns:
        sample_id, HLA_locus, HLA_allele, vdjdb_antigen,
        vdjdb_antigen_species, n_clonotypes, n_cells

Usage:
    python summarize_hla_epitopes.py \
        --input sample_clonotype_summary.vdjdb.tsv \
        --output sample_hla_epitope_summary.tsv
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from typing import Dict, List, Tuple

import pandas as pd

HLA_LOCI: List[str] = ["A", "B", "C", "DRB1", "DQB1", "DQA1", "DPA1", "DPB1"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize annotated clonotypes by HLA allele and antigen."
    )
    parser.add_argument("--input", required=True, help="Annotated clonotype TSV.")
    parser.add_argument("--output", required=True, help="Output summary TSV.")
    return parser.parse_args()


def _extract_vdjdb_entries(row: pd.Series) -> List[Tuple[str, str]]:
    """
    Returns list of (antigen, species) tuples from vdjdb_entries JSON column.
    Falls back to splitting vdjdb_antigen/species strings if JSON absent.
    """
    entries_json = row.get("vdjdb_entries", "")
    if isinstance(entries_json, str) and entries_json.strip().startswith("["):
        try:
            entries = json.loads(entries_json)
            return [(e[0], e[1]) for e in entries]
        except json.JSONDecodeError:
            pass
    antigens = str(row.get("vdjdb_antigen", "") or "")
    species = str(row.get("vdjdb_antigen_species", "") or "")
    antigen_list = [a for a in antigens.split(";") if a]
    species_list = [s for s in species.split(";") if s]
    if not antigen_list:
        return []
    if len(species_list) != len(antigen_list):
        species_list = [""] * len(antigen_list)
    return list(zip(antigen_list, species_list))


def summarize(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path, sep="\t")
    required = ["sample_id", "vdjdb_match", "n_cells"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Annotated clonotype summary missing columns: {missing}")
    df = df[df["vdjdb_match"].astype(str).str.upper().isin(["TRUE", "1", "YES"])]
    if df.empty:
        pd.DataFrame(
            columns=[
                "sample_id",
                "HLA_locus",
                "HLA_allele",
                "vdjdb_antigen",
                "vdjdb_antigen_species",
                "n_clonotypes",
                "n_cells",
            ]
        ).to_csv(output_path, sep="\t", index=False)
        return

    summary_counter: Dict[
        Tuple[str, str, str, str, str], Dict[str, float]
    ] = defaultdict(lambda: {"n_clonotypes": 0, "n_cells": 0.0})

    for _, row in df.iterrows():
        entries = _extract_vdjdb_entries(row)
        if not entries:
            continue
        n_cells_raw = row.get("n_cells", 0)
        n_cells = 0.0 if pd.isna(n_cells_raw) else float(n_cells_raw)
        used_clonotype = set()
        for locus in HLA_LOCI:
            col1 = f"HLA_{locus}_allele1"
            col2 = f"HLA_{locus}_allele2"
            for allele_col in (col1, col2):
                allele = row.get(allele_col)
                if pd.isna(allele) or not allele:
                    continue
                for antigen, species in entries:
                    key = (
                        row["sample_id"],
                        locus,
                        str(allele),
                        antigen,
                        species,
                    )
                    if key not in used_clonotype:
                        summary_counter[key]["n_clonotypes"] += 1
                        summary_counter[key]["n_cells"] += n_cells
                        used_clonotype.add(key)

    records: List[Dict[str, object]] = []
    for key, metrics in summary_counter.items():
        sample_id, locus, allele, antigen, species = key
        records.append(
            {
                "sample_id": sample_id,
                "HLA_locus": locus,
                "HLA_allele": allele,
                "vdjdb_antigen": antigen,
                "vdjdb_antigen_species": species,
                "n_clonotypes": int(metrics["n_clonotypes"]),
                "n_cells": metrics["n_cells"],
            }
        )
    summary_df = pd.DataFrame(records)
    summary_df.sort_values(
        ["sample_id", "HLA_locus", "vdjdb_antigen", "HLA_allele"], inplace=True
    )
    summary_df.to_csv(output_path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    summarize(args.input, args.output)


if __name__ == "__main__":
    main()

