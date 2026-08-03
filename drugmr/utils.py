#!/usr/bin/env python3
# utils.py
import polars as pl


def filter_mr_targets(df: pl.DataFrame):
    targets = []
    for row in df.iter_rows(named=True):
        protein = row["protein"]
        n_instruments = row["n_instruments"]
        wald = row["Wald_FDR_q"]
        q = row["Q_pval"]
        ivw = row["IVW_FDR_q"]
        if n_instruments == 1:
            if wald is not None and wald < 0.05:
                targets.append(protein)
        elif n_instruments > 1:
            if (ivw is not None and q is not None and ivw < 0.05 and q > 0.05):
                targets.append(protein)
    return targets