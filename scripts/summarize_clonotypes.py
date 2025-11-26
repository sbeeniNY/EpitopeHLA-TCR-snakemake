#!/usr/bin/env python3
"""
Clonotype-level summarization.

This script ingests a merged TCR-HLA TSV produced by merge_tcr_hla.py and
builds a clonotype-level summary table. It groups rows by sample_id,
chain, and exact_subclonotype_id (falling back to raw_clonotype_id when
the exact ID is missing), aggregates key statistics, and carries over
HLA genotype columns for QC.

Expected input columns (tab-delimited UTF-8):
    sample_id, barcode, chain, exact_subclonotype_id, raw_clonotype_id,
    cdr3, cdr3_nt, v_gene, j_gene, reads, umis,
    HLA_*_allele{1,2} (e.g. HLA_A_allele1).

Output TSV columns (tab-delimited UTF-8):
    sample_id, clonotype_id, chain, n_cells, n_contigs,
    cdr3_aa, cdr3_nt, v_gene, j_gene, total_reads, total_umis,
    HLA_*_allele1/2 ...

This script is designed to be called from Snakemake, but it can be run
manually:
    python summarize_clonotypes.py --input sample_tcr_hla.tsv \
        --output sample_clonotype_summary.tsv
"""
from __future__ import annotations

import argparse
import logging
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

HLA_LOCI: List[str] = ["A", "B", "C", "DRB1", "DQB1", "DQA1", "DPA1", "DPB1"]
REPR_TIEBREAK_COLUMNS = ["reads", "umis"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize merged TCR+HLA TSV into clonotype-level table."
    )
    parser.add_argument("--input", required=True, help="Merged TCR-HLA TSV path.")
    parser.add_argument("--output", required=True, help="Output clonotype summary TSV.")
    return parser.parse_args()


def validate_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input file missing required columns: {missing}")


def determine_clonotype_id(df: pd.DataFrame) -> pd.Series:
    has_exact = "exact_subclonotype_id" in df.columns
    has_raw = "raw_clonotype_id" in df.columns
    if not has_exact and not has_raw:
        raise ValueError(
            "Input must contain either 'exact_subclonotype_id' or 'raw_clonotype_id'."
        )
    exact = df["exact_subclonotype_id"] if has_exact else pd.Series([pd.NA] * len(df))
    raw = df["raw_clonotype_id"] if has_raw else pd.Series([pd.NA] * len(df))
    clonotype = exact.fillna(raw)
    if clonotype.isna().any():
        raise ValueError(
            "Clonotype identifier missing even after fallback to raw_clonotype_id."
        )
    return clonotype.astype(str)


def _pick_representative(
    subdf: pd.DataFrame, column: str, tiebreak_cols: Iterable[str]
) -> Optional[str]:
    if column not in subdf.columns:
        return None
    values = subdf[column].dropna()
    if values.empty:
        return None
    counts = Counter(values)
    top_freq = max(counts.values())
    candidates = [val for val, freq in counts.items() if freq == top_freq]
    if len(candidates) == 1:
        return candidates[0]
    # Tie-break using summed weights from specified columns.
    best_val = candidates[0]
    best_weight = float("-inf")
    for val in candidates:
        mask = subdf[column] == val
        weight = 0.0
        for wcol in tiebreak_cols:
            if wcol in subdf.columns:
                weight += subdf.loc[mask, wcol].fillna(0).astype(float).sum()
        if weight > best_weight:
            best_weight = weight
            best_val = val
    return best_val


def _summarize_group(subdf: pd.DataFrame) -> Dict[str, object]:
    result: Dict[str, object] = {}
    barcode_col = "barcode" if "barcode" in subdf.columns else None
    result["n_cells"] = (
        subdf[barcode_col].nunique(dropna=True) if barcode_col else pd.NA
    )
    result["n_contigs"] = len(subdf)
    reads_col = "reads" if "reads" in subdf.columns else None
    umis_col = "umis" if "umis" in subdf.columns else None
    result["total_reads"] = (
        subdf[reads_col].fillna(0).astype(float).sum() if reads_col else 0
    )
    result["total_umis"] = (
        subdf[umis_col].fillna(0).astype(float).sum() if umis_col else 0
    )
    result["cdr3_aa"] = _pick_representative(
        subdf, "cdr3", REPR_TIEBREAK_COLUMNS
    ) or _pick_representative(subdf, "cdr3_aa", REPR_TIEBREAK_COLUMNS)
    result["cdr3_nt"] = _pick_representative(subdf, "cdr3_nt", REPR_TIEBREAK_COLUMNS)
    result["v_gene"] = _pick_representative(subdf, "v_gene", REPR_TIEBREAK_COLUMNS)
    result["j_gene"] = _pick_representative(subdf, "j_gene", REPR_TIEBREAK_COLUMNS)
    return result


def _collect_hla_fields(subdf: pd.DataFrame) -> Dict[str, object]:
    hla_values: Dict[str, object] = {}
    for locus in HLA_LOCI:
        col1 = f"HLA_{locus}_allele1"
        col2 = f"HLA_{locus}_allele2"
        if col1 not in subdf.columns or col2 not in subdf.columns:
            continue
        combos = (
            subdf[[col1, col2]]
            .dropna(how="all")
            .astype(str)
            .drop_duplicates()
        )
        if len(combos) == 0:
            hla_values[col1] = pd.NA
            hla_values[col2] = pd.NA
        elif len(combos) == 1:
            hla_values[col1] = combos.iloc[0, 0]
            hla_values[col2] = combos.iloc[0, 1]
        else:
            logging.warning(
                "Multiple HLA allele combinations detected for locus %s: %s. "
                "Using most frequent combination.",
                locus,
                combos.values.tolist(),
            )
            combo_counts = (
                subdf[[col1, col2]]
                .dropna(how="all")
                .astype(str)
                .value_counts()
            )
            top = combo_counts.idxmax()
            hla_values[col1] = top[0]
            hla_values[col2] = top[1]
    return hla_values


def summarize(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path, sep="\t")
    validate_columns(df, ["sample_id", "chain"])
    df["clonotype_id"] = determine_clonotype_id(df)
    group_cols = ["sample_id", "chain", "clonotype_id"]
    rows: List[Dict[str, object]] = []
    for (sample_id, chain, clonotype_id), subdf in df.groupby(group_cols, dropna=False):
        summary = _summarize_group(subdf)
        summary.update(
            {
                "sample_id": sample_id,
                "chain": chain,
                "clonotype_id": clonotype_id,
            }
        )
        summary.update(_collect_hla_fields(subdf))
        rows.append(summary)
    result_df = pd.DataFrame(rows)
    # Ensure consistent column ordering.
    base_cols = [
        "sample_id",
        "clonotype_id",
        "chain",
        "n_cells",
        "n_contigs",
        "cdr3_aa",
        "cdr3_nt",
        "v_gene",
        "j_gene",
        "total_reads",
        "total_umis",
    ]
    hla_cols: List[str] = []
    for locus in HLA_LOCI:
        col1 = f"HLA_{locus}_allele1"
        col2 = f"HLA_{locus}_allele2"
        if col1 in result_df.columns:
            hla_cols.append(col1)
        if col2 in result_df.columns:
            hla_cols.append(col2)
    ordered_cols = base_cols + hla_cols
    for col in ordered_cols:
        if col not in result_df.columns:
            result_df[col] = pd.NA
    result_df = result_df[ordered_cols]
    result_df.sort_values(["sample_id", "chain", "clonotype_id"], inplace=True)
    result_df.to_csv(output_path, sep="\t", index=False)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    summarize(args.input, args.output)


if __name__ == "__main__":
    main()

