#!/usr/bin/env python3
"""
Annotate clonotype summaries with VDJdb information.

Inputs:
    - Clonotype summary TSV from summarize_clonotypes.py
    - VDJdb-like TSV (e.g. resources/vdjdb/vdjdb.txt)

Outputs:
    - Clonotype summary TSV with additional VDJdb annotation columns:
        vdjdb_match (TRUE/FALSE)
        vdjdb_n_matches
        vdjdb_antigen
        vdjdb_antigen_species
        vdjdb_mhc
        vdjdb_reference

Matching logic:
    - Exact match on chain (e.g. TRB) and CDR3 amino acid sequence.
    - Aggregates all matches (unique tuples) and joins values with ';'.
    - Column mapping from VDJdb -> script is configurable via CLI.

Example:
    python annotate_vdjdb.py --input sample_clonotype_summary.tsv \
        --vdjdb resources/vdjdb/vdjdb.txt \
        --output sample_clonotype_summary.vdjdb.tsv
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import pandas as pd


@dataclass(frozen=True)
class VDJDBColumns:
    cdr3: str = "cdr3"
    gene: str = "gene"
    v_gene: str = "v"
    j_gene: str = "j"
    antigen: str = "antigen.epitope"
    antigen_species: str = "antigen.species"
    mhc: str = "mhc.a"
    reference: str = "reference.id"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate clonotype summary with VDJdb matches."
    )
    parser.add_argument("--input", required=True, help="Clonotype summary TSV.")
    parser.add_argument("--vdjdb", required=True, help="VDJdb TSV file.")
    parser.add_argument("--output", required=True, help="Annotated TSV output path.")
    parser.add_argument(
        "--vdjdb-cdr3-col",
        default="cdr3",
        help="Column in VDJdb file for CDR3 amino acid sequence (default: cdr3).",
    )
    parser.add_argument(
        "--vdjdb-gene-col",
        default="gene",
        help="Column in VDJdb file for chain/gene (default: gene).",
    )
    parser.add_argument(
        "--vdjdb-antigen-col",
        default="antigen.epitope",
        help="Column in VDJdb file for antigen/epitope (default: antigen.epitope).",
    )
    parser.add_argument(
        "--vdjdb-antigen-species-col",
        default="antigen.species",
        help="Column for antigen species (default: antigen.species).",
    )
    parser.add_argument(
        "--vdjdb-mhc-col",
        default="mhc.a",
        help="Column for MHC/HLA restriction (default: mhc.a).",
    )
    parser.add_argument(
        "--vdjdb-reference-col",
        default="reference.id",
        help="Column for reference/PMID (default: reference.id).",
    )
    return parser.parse_args()


def _prepare_vdjdb(
    vdjdb_path: str,
    columns: VDJDBColumns,
) -> Dict[Tuple[str, str], Dict[str, List[str]]]:
    df = pd.read_csv(vdjdb_path, sep="\t", dtype=str)
    required = [columns.cdr3, columns.gene]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"VDJdb file missing required columns: {missing}")
    df = df.copy()
    df[columns.cdr3] = df[columns.cdr3].str.upper().str.strip()
    df[columns.gene] = df[columns.gene].str.upper().str.strip()
    df.dropna(subset=[columns.cdr3, columns.gene], inplace=True)
    annotation_map: Dict[Tuple[str, str], Dict[str, List[str]]] = {}
    for (gene, cdr3), subdf in df.groupby([columns.gene, columns.cdr3]):
        entries = []
        for _, row in subdf.iterrows():
            entries.append(
                {
                    "antigen": str(row.get(columns.antigen, "") or ""),
                    "antigen_species": str(row.get(columns.antigen_species, "") or ""),
                    "mhc": str(row.get(columns.mhc, "") or ""),
                    "reference": str(row.get(columns.reference, "") or ""),
                }
            )
        annotation_map[(gene, cdr3)] = {
            "entries": entries,
            "count": len(entries),
        }
    return annotation_map


def _aggregate_entries(entries: List[Dict[str, str]]) -> Dict[str, str]:
    # Use tuple to enforce consistent ordering when joining.
    unique_entries = []
    seen = set()
    for entry in entries:
        key = (
            entry.get("antigen", ""),
            entry.get("antigen_species", ""),
            entry.get("mhc", ""),
            entry.get("reference", ""),
        )
        if key not in seen:
            seen.add(key)
            unique_entries.append(key)
    unique_entries.sort()
    antigens = ";".join(item[0] for item in unique_entries if item[0])
    species = ";".join(item[1] for item in unique_entries if item[1])
    mhc = ";".join(item[2] for item in unique_entries if item[2])
    references = ";".join(item[3] for item in unique_entries if item[3])
    return {
        "vdjdb_antigen": antigens,
        "vdjdb_antigen_species": species,
        "vdjdb_mhc": mhc,
        "vdjdb_reference": references,
        "vdjdb_entries": json.dumps(unique_entries, ensure_ascii=False),
    }


def annotate_clonotypes(
    clonotype_path: str,
    vdjdb_map: Dict[Tuple[str, str], Dict[str, List[str]]],
    output_path: str,
) -> None:
    df = pd.read_csv(clonotype_path, sep="\t")
    required = ["chain", "cdr3_aa", "n_cells"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Clonotype summary missing required columns: {missing}")

    match_cols = {
        "vdjdb_match": [],
        "vdjdb_n_matches": [],
        "vdjdb_antigen": [],
        "vdjdb_antigen_species": [],
        "vdjdb_mhc": [],
        "vdjdb_reference": [],
        "vdjdb_entries": [],
    }

    for _, row in df.iterrows():
        chain = str(row["chain"]).upper().strip()
        cdr3 = str(row["cdr3_aa"]).upper().strip()
        key = (chain, cdr3)
        if cdr3 in ("", "NAN") or key not in vdjdb_map:
            match_cols["vdjdb_match"].append(False)
            match_cols["vdjdb_n_matches"].append(0)
            match_cols["vdjdb_antigen"].append("")
            match_cols["vdjdb_antigen_species"].append("")
            match_cols["vdjdb_mhc"].append("")
            match_cols["vdjdb_reference"].append("")
            match_cols["vdjdb_entries"].append("[]")
        else:
            info = vdjdb_map[key]
            aggregated = _aggregate_entries(info["entries"])
            match_cols["vdjdb_match"].append(True)
            match_cols["vdjdb_n_matches"].append(info["count"])
            match_cols["vdjdb_antigen"].append(aggregated["vdjdb_antigen"])
            match_cols["vdjdb_antigen_species"].append(
                aggregated["vdjdb_antigen_species"]
            )
            match_cols["vdjdb_mhc"].append(aggregated["vdjdb_mhc"])
            match_cols["vdjdb_reference"].append(aggregated["vdjdb_reference"])
            match_cols["vdjdb_entries"].append(aggregated["vdjdb_entries"])

    for col, values in match_cols.items():
        df[col] = values

    df.to_csv(output_path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    columns = VDJDBColumns(
        cdr3=args.vdjdb_cdr3_col,
        gene=args.vdjdb_gene_col,
        antigen=args.vdjdb_antigen_col,
        antigen_species=args.vdjdb_antigen_species_col,
        mhc=args.vdjdb_mhc_col,
        reference=args.vdjdb_reference_col,
    )
    vdjdb_map = _prepare_vdjdb(args.vdjdb, columns)
    annotate_clonotypes(args.input, vdjdb_map, args.output)


if __name__ == "__main__":
    main()

