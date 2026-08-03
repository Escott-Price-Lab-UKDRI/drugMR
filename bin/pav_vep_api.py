#!/usr/bin/env python3
import polars as pl
import requests
import argparse
from pathlib import Path
from drugmr.utils import filter_mr_targets


# prospective plan 
# for each cis-MR significant target in pQTL dataset X 
# map their genetic instruments from instruments/ and also A1 and A2
# annotate those and retain in temp_file within work/ 
# for each significant target - find cis-region with respect to pQTL dataset 
# annotate all of them and where for SNP X -> GENE_id > 1 -> retain where A1 & A1 == cannonical
# for target X (not save all files) -> do it per protein 
# map all snps in pav_terms - check whether either == pav term - if so new col in work/ dataset -> boolean col
# for each SNP which == PAV within target X -> r2 LD (with PLINK) of instruments vs that SNP
# DONE FOR NOW


pav_terms = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

def vep_api_request(targets: dict[str, str]) -> pl.DataFrame:
    r = requests.post(
        "https://rest.ensembl.org/vep/human/id",
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json={"ids": list(targets), "mane": 1, "hgvs": 1, "protein": 1},
    )
    r.raise_for_status()
    rows = [
        {
            "variant_id": v["input"],
            "variant_allele": tc.get("variant_allele"),
            "gene_symbol": tc.get("gene_symbol"),
            "gene_id": tc.get("gene_id"),
            "transcript_id": tc.get("transcript_id"),
            "consequence": tc.get("consequence_terms", []),
            "hgvsc": tc.get("hgvsc"),
            "hgvsp": tc.get("hgvsp"),
        }
        for v in r.json()
        for tc in v.get("transcript_consequences", [])
        if tc.get("mane_select")
        and tc.get("gene_symbol") == targets.get(v["input"])
    ]

    if not rows:
        return pl.DataFrame()
    return (pl.DataFrame(rows).explode("consequence", empty_as_null=True).with_columns(pl.col("consequence").is_in(pav_terms).alias("is_pav")))


def pav_vep_checks(pqtl_dataset: str, pheno_id: str):
    # ukb_ppp_AD_all_MR_instruments.tsv
    # wingo_brain_AD_all_MR.tsv
    cis_mr_instruments = Path(f"./results/cis-MR/instruments/{pqtl_dataset}_{pheno_id}_all_MR_instruments.tsv")
    cis_mr_res = Path(f"./results/cis-MR/{pqtl_dataset}_{pheno_id}_all_MR.tsv")
    cis_mr_instruments = pl.read_csv(cis_mr_instruments, separator="\t")
    cis_mr_res = pl.read_csv(cis_mr_res, separator="\t")
    # filter significant MR targets
    candidates = filter_mr_targets(cis_mr_res)
    processed = set()
    vep_results = []

    print(f"MR candidates: {len(candidates)}")
    for row in cis_mr_instruments.iter_rows(named=True):
        protein = row["protein"]
        if protein in candidates and protein not in processed:
            processed.add(protein)
            # all SNPs used as instruments for that target
            temp_df = cis_mr_instruments.filter(pl.col("protein") == protein)
            gene = protein.split("_")[0]
            ivs_to_targets = {}

            for snp in temp_df["SNP"]:
                ivs_to_targets[snp] = gene

            print(f"Checking {protein}: {len(ivs_to_targets)} instrument(s)")

            if ivs_to_targets:
                result = vep_api_request(targets=ivs_to_targets)

                if not result.is_empty():
                    print(result)
                    vep_results.append(result)
                else:
                    print(f"No matching MANE consequence for {protein}")

    if vep_results:
        vep_results = pl.concat(vep_results, how="diagonal_relaxed")
    else:
        vep_results = pl.DataFrame()

    print("\nCombined VEP results:")
    print(vep_results)
    return vep_results


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno-id", "--pheno_id", dest="pheno_id", required=True)
    p.add_argument("--pqtl-dataset", "--pqtl_dataset", dest="pqtl_dataset", required=True)
    args = p.parse_args()
    pav_vep_checks(
        pheno_id=args.pheno_id,
        pqtl_dataset=args.pqtl_dataset,
    )


if __name__ == "__main__":
    main()