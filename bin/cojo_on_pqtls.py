#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path
import os
from drugmr import COJO

def cojo_on_pqtl_loci(ref_bfile: str, pqtl_dataset: str, pheno_id: str):
    
    """
    1. Check proteins within pairwise pQTL-GWAS COLOC results for pQTL dataset X
    2. Reformat pqtl.parquet in cis-regions/ to COJO and store as temp in work/COJO
    3. Align all pQTL effects to outcome GWAS A1
    4. Run COJO
    """

    # all proteins in this file have already passed cis-MR + COLOC thresholds
    res_file = f"./results/coloc/{pqtl_dataset}/{pqtl_dataset}_{pheno_id}_all_coloc.tsv"
    df = pl.read_csv(res_file, separator="\t")
    # coloc.R names the protein column protein_id
    if "protein_id" in df.columns:
        df = df.rename({"protein_id": "protein"})
    print(f"[TRACKING] Proteins available for COJO: {df.height}")
    temp_out_path = f"./work/COJO/{pqtl_dataset}/{pheno_id}"
    temp_out_path = Path(temp_out_path)
    os.makedirs(temp_out_path, exist_ok=True)
    # bare in mind top SNP
    for row in df.iter_rows(named=True):
        target = row["protein"]
        pth = f"./dat/cis_regions/{pqtl_dataset}/{target}"
        pth = Path(pth)
        target_cis_region = pth / "pqtl.parquet"
        target_gwas_region = pth / "gwas.parquet"
        # print(f"\n[{i}/{df.height}] {target}")
        print(f"Reading: {target_cis_region}")
        cis = pl.read_parquet(target_cis_region)
        print(f"Loaded {cis.height:,} pQTL SNPs")
        print("pQTL columns:", cis.columns)
        # read outcome GWAS cis-region
        # GWAS A1 == effect allele used across all downstream analyses
        print(f"Reading: {target_gwas_region}")
        gwas = pl.read_parquet(target_gwas_region)
        print(f"Loaded {gwas.height:,} GWAS SNPs")
        print("GWAS columns:", gwas.columns)
        # standardise allele columns
        cis = cis.with_columns(pl.col("SNP").cast(pl.Utf8), pl.col("A1").cast(pl.Utf8).str.to_uppercase(), pl.col("A2").cast(pl.Utf8).str.to_uppercase())
        gwas = (gwas.select(pl.col("SNP").cast(pl.Utf8), pl.col("A1").cast(pl.Utf8).str.to_uppercase().alias("GWAS_A1"), pl.col("A2").cast(pl.Utf8).str.to_uppercase().alias("GWAS_A2"), pl.col("BETA").alias("GWAS_BETA")).unique(subset=["SNP"], keep="first"))
        # match pQTL SNPs against outcome GWAS SNPs
        cis = cis.join(gwas, on="SNP", how="inner")
        print(f"Retained {cis.height:,} SNPs present in both pQTL and GWAS files")
        # align pQTL BETA + FRQ to outcome GWAS A1
        cis = cis.with_columns(
            pl.when(
                (pl.col("A1") == pl.col("GWAS_A1")) &
                (pl.col("A2") == pl.col("GWAS_A2"))
            )
            .then(pl.col("BETA"))
            .when(
                (pl.col("A1") == pl.col("GWAS_A2")) &
                (pl.col("A2") == pl.col("GWAS_A1"))
            )
            .then(-pl.col("BETA"))
            .otherwise(None)
            .alias("BETA_ALIGNED"),

            pl.when(
                (pl.col("A1") == pl.col("GWAS_A1")) &
                (pl.col("A2") == pl.col("GWAS_A2"))
            )
            .then(pl.col("FRQ"))
            .when(
                (pl.col("A1") == pl.col("GWAS_A2")) &
                (pl.col("A2") == pl.col("GWAS_A1"))
            )
            .then(1 - pl.col("FRQ"))
            .otherwise(None)
            .alias("A1FREQ_ALIGNED"),

            pl.when(
                (pl.col("A1") == pl.col("GWAS_A1")) &
                (pl.col("A2") == pl.col("GWAS_A2"))
            )
            .then(pl.lit(False))
            .when(
                (pl.col("A1") == pl.col("GWAS_A2")) &
                (pl.col("A2") == pl.col("GWAS_A1"))
            )
            .then(pl.lit(True))
            .otherwise(None)
            .alias("PQTl_FLIPPED"),
        )

        n_before_alignment = cis.height
        cis = cis.drop_nulls(["GWAS_A1", "GWAS_A2", "BETA_ALIGNED", "A1FREQ_ALIGNED"])
        print(
            f"Retained {cis.height:,}/{n_before_alignment:,} SNPs "
            f"after aligning pQTL effects to {pheno_id} GWAS A1"
        )

        # reformat cis-pQTL summary statistics for COJO
        # ALLELE1 == outcome GWAS A1
        # ALLELE0 == outcome GWAS A2
        # BETA == pQTL effect aligned to outcome GWAS A1
        cis_cojo = cis.select(
            pl.col("SNP").alias("ID"),
            pl.col("GWAS_A1").alias("ALLELE1"),
            pl.col("GWAS_A2").alias("ALLELE0"),
            pl.col("A1FREQ_ALIGNED").alias("A1FREQ"),
            pl.col("BETA_ALIGNED").alias("BETA"),
            pl.col("SE"),
            pl.col("P"),
            pl.col("N"),
        )

        # check the lead SNP alignment
        lead_snp = row["top_snp"]
        lead_snp_check = (
            cis
            .filter(pl.col("SNP") == lead_snp)
            .select(
                "SNP",
                "A1",
                "A2",
                "GWAS_A1",
                "GWAS_A2",
                "GWAS_BETA",
                "BETA",
                "BETA_ALIGNED",
                "FRQ",
                "A1FREQ_ALIGNED",
                "PQTl_FLIPPED",
            )
        )

        if lead_snp_check.height > 0:
            print(f"Lead SNP alignment for {target}:")
            print(lead_snp_check)
        else:
            print(
                f"Lead SNP {lead_snp} was not retained after "
                f"pQTL-GWAS allele alignment for {target}"
            )

        # target_cojo.input
        cojo_file = temp_out_path / f"{target}_cojo.input"
        cis_cojo.write_csv(cojo_file, separator="\t")
        print(f"Wrote {cojo_file}")
        # COJO output prefix
        out_prefix = f"./results/COJO/{pqtl_dataset}/{pheno_id}/{target}/{target}"
        out_prefix = Path(out_prefix)
        os.makedirs(out_prefix.parent, exist_ok=True)
        # run COJO
        print("Running COJO...")
        COJO(
            cojo_sumstats=cojo_file,
            p_thresh=5e-8,
            ref_bfile=ref_bfile,
            collinear_thresh=0.9,
            wind_thresh=10000,
            out_prefix=out_prefix,
        )
        # results/COJO/ukb_ppp/AD/BLNK_Q8WV28/BLNK_Q8WV28.jma.cojo

        # now check COJO output for independent signals
        cojo_max_signals = 10
        mr_instruments_file = Path(f"./results/cis-MR/instruments/{pqtl_dataset}_{pheno_id}_all_MR_instruments.tsv")
        mr_instruments = pl.read_csv(mr_instruments_file, separator="\t", schema_overrides={"protein": pl.Utf8, "SNP": pl.Utf8})
        cojo_results_file = Path(f"{out_prefix}.jma.cojo")
        cojo_results = pl.read_csv(cojo_results_file, separator="\t", schema_overrides={"SNP": pl.Utf8})
        n_signals = cojo_results.height
        print(f"[TRACKING] COJO selected {n_signals} independent signal(s) for {target}")
        if n_signals > cojo_max_signals:
            print(f"[CONCERN] {target} has {n_signals} COJO signals (> {cojo_max_signals}). Keeping only SNPs used in cis-MR.")
            target_mr_snps = mr_instruments.filter(pl.col("protein") == target).select("SNP").unique()
            cojo_results = cojo_results.join(target_mr_snps, on="SNP", how="semi")
            cojo_results.write_csv(cojo_results_file, separator="\t")
            print(f"[TRACKING] Overwrote {cojo_results_file} with {cojo_results.height} COJO rows matching cis-MR instruments")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--ref_bfile", required=True)
    args = p.parse_args()
    cojo_on_pqtl_loci(
        ref_bfile=args.ref_bfile,
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
    )

if __name__ == "__main__":
    main()