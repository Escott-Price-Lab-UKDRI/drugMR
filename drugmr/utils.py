#!/usr/bin/env python3
# utils.py
import polars as pl
import subprocess


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

# plink \
#   --bfile 1000G.EUR.QC.ALL \
#   --ld rs653765 rs1427281

def impute_ld(ref_bfile, snp_1, snp_2):
    cmd = f"""
plink \
    --bfile {ref_bfile} \
    --ld {snp_1} {snp_2}
"""
    return subprocess.run(cmd, shell=True, check=False, executable="/bin/bash", capture_output=True, text=True)