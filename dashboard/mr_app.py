#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st


# KEY CHANGES DOWN THE LINE WITH MORE PQTL DATASETS 
# -> CHANGE THE DASHBOARD FUNCT TO ADD MORE PQTL DATASETS
# biomarker meta analysis: https://pmc.ncbi.nlm.nih.gov/articles/instance/12136742/pdf/nihpp-rs6597595v1.pdf

def create_streamlit_ammenities(db_name: str, port_number: str):
    cmd = f"""
set -euo pipefail 
mkdir -p .streamlit
cat > .streamlit/secrets.toml <<EOF
[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "{port_number}"
database = "{db_name}"
username = ""
password = ""
EOF
    """

    # run in terminal to create streamlit ammenities 
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def retention(current: int, previous: int):
    return 0.0 if previous == 0 else 100 * current / previous


def available_cols(df: pd.DataFrame, cols: list[str]):
    return [col for col in cols if col in df.columns]


def load_required_tsv(file: Path, label: str):
    if not file.exists():
        st.error(f"{label} result file not found: {file}")
        st.stop()

    df = pd.read_csv(file, sep="\t")

    if df.empty:
        st.error(f"{label} result file is empty: {file}")
        st.stop()

    return df


def load_optional_tsv(file: Path, label: str):
    if not file.exists():
        st.warning(f"{label} result file not found: {file}")
        return pd.DataFrame()

    df = pd.read_csv(file, sep="\t")

    if df.empty:
        st.warning(f"{label} result file is empty: {file}")
        return pd.DataFrame()

    return df


def find_result_file(project_dir: Path, candidate_files: list[Path], candidate_names: list[str]):
    for file in candidate_files:
        if file.exists():
            return file

    matches = []

    for candidate_name in candidate_names:
        matches.extend(list((project_dir / "results").rglob(candidate_name)))

    matches = sorted(set(matches))

    if len(matches) == 1:
        return matches[0]

    return None


def filter_protein(df: pd.DataFrame, protein: str):
    if df.empty or not protein or "protein" not in df.columns:
        return df

    return df[
        df["protein"]
        .astype(str)
        .str.contains(protein, case=False, na=False)
    ].copy()


def safe_nunique(df: pd.DataFrame, col: str):
    if df.empty or col not in df.columns:
        return 0

    return df[col].nunique()


def standardise_columns(df: pd.DataFrame):
    df = df.copy()

    df.columns = (
        pd.Index(df.columns)
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"[^a-z0-9]+", "_", regex=True)
        .str.strip("_")
    )

    return df



def dashboard(db_name: str, port_number: str, phenotype: str, pqtl_dataset: str):
    mr_table = "cis_mr_results"
    coloc_table = "coloc_results"
    phewas_table = "phewas_safety"

    # main aesthetics
    st.set_page_config(page_title=f"{db_name}", layout="wide")
    st.markdown(
        """
        <style>
        .block-container {padding-top: 2rem; padding-bottom: 3rem;}
        [data-testid="stMetric"] {border: 1px solid rgba(49, 51, 63, 0.15); border-radius: 10px; padding: 14px 16px;}
        [data-testid="stDataFrame"] {border: 1px solid rgba(49, 51, 63, 0.12); border-radius: 10px; overflow: hidden;}
        div[data-baseweb="select"] > div {border-radius: 8px;}
        </style>
        """,
        unsafe_allow_html=True
    )
    conn = st.connection("postgresql", type="sql", url=f"postgresql://localhost:{port_number}/{db_name}")

    # pQTL dataset selection schema 
    # CLI pQTL dataset is used as the default dashboard selection
    dataset_names = {
        "ukb_ppp": "UKB-PPP",
        "decode": "deCODE",
        "wu_csf": "WU-CSF"
    }

    dataset_ns = {
        "ukb_ppp": 54219,
        "decode": 35559,
        "wu_csf": 3506
    }

    project_dir = Path(__file__).resolve().parent.parent

    # check which datasets have the required dashboard files
    dataset_result_files = {}
    available_datasets = []

    for dataset_id in dataset_names:
        mr_file = find_result_file(
            project_dir,
            [
                project_dir / f"results/cis-MR/{dataset_id}_{phenotype}_all_MR.tsv",
                project_dir / f"results/cis_MR/{dataset_id}_{phenotype}_all_MR.tsv",
                project_dir / f"results/cis_MR/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_MR.tsv",
                project_dir / f"results/MR/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_MR.tsv",
                project_dir / f"results/MR/{dataset_id}_{phenotype}_all_MR.tsv"
            ],
            [
                f"{dataset_id}_{phenotype}_all_MR.tsv",
                f"{dataset_id}_{phenotype}_MR.tsv"
            ]
        )

        coloc_file = find_result_file(
            project_dir,
            [
                project_dir / f"results/coloc/{dataset_id}/{dataset_id}_{phenotype}_all_coloc.tsv",
                project_dir / f"results/coloc/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_coloc.tsv",
                project_dir / f"results/COLOC/{dataset_id}/{phenotype}/{dataset_id}_{phenotype}_all_coloc.tsv",
                project_dir / f"results/COLOC/{dataset_id}_{phenotype}_all_coloc.tsv",
                project_dir / f"results/coloc/{dataset_id}_{phenotype}_all_coloc.tsv"
            ],
            [
                f"{dataset_id}_{phenotype}_all_coloc.tsv",
                f"{dataset_id}_{phenotype}_coloc.tsv",
                f"{dataset_id}_{phenotype}_COLOC.tsv"
            ]
        )

        phewas_file = project_dir / "results" / "PheWAS" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_PheWAS.tsv"
        target_info_file = project_dir / "results" / "target_stats" / dataset_id / phenotype / f"{dataset_id}_{phenotype}_top_cis_hits.tsv"

        required_files = [
            mr_file,
            coloc_file
        ]

        if all(file is not None and file.exists() for file in required_files):
            available_datasets.append(dataset_id)
            dataset_result_files[dataset_id] = {
                "mr": mr_file,
                "coloc": coloc_file,
                "phewas": phewas_file,
                "target_info": target_info_file
            }

    if len(available_datasets) == 0:
        st.error(f"No dataset has a complete set of cis-MR and COLOC dashboard files for {phenotype}.")
        st.stop()

    # use the CLI dataset as default
    # otherwise use the first complete dataset which was found
    if pqtl_dataset not in available_datasets:
        pqtl_dataset = available_datasets[0]

    available_dataset_names = [dataset_names[dataset_id] for dataset_id in available_datasets]
    default_dataset_name = dataset_names[pqtl_dataset]
    selected_dataset_name = st.segmented_control(
        "pQTL dataset",
        available_dataset_names,
        default=default_dataset_name,
        selection_mode="single",
        key="pqtl_dataset_selector"
    )

    dataset_ids = {dataset_name: dataset_id for dataset_id, dataset_name in dataset_names.items()}
    pqtl_dataset = dataset_ids[selected_dataset_name]
    dataset_name = dataset_names[pqtl_dataset]
    dataset_n = dataset_ns[pqtl_dataset]

    # corresponding selected dataset result files
    mr_file = dataset_result_files[pqtl_dataset]["mr"]
    coloc_file = dataset_result_files[pqtl_dataset]["coloc"]
    phewas_file = dataset_result_files[pqtl_dataset]["phewas"]
    target_info_file = dataset_result_files[pqtl_dataset]["target_info"]

    # load local result files into PostgreSQL for the dashboard
    mr = load_required_tsv(mr_file, "cis-MR")
    coloc = load_required_tsv(coloc_file, "pQTL–GWAS COLOC")
    phewas = load_optional_tsv(phewas_file, "PheWAS safety")
    target_info = load_optional_tsv(target_info_file, "Harmonised target information")

    # standardise MR + pQTL COLOC columns before loading into PostgreSQL
    # avoids dataset-specific differences such as Wald_beta vs wald_beta
    mr = standardise_columns(mr)
    coloc = standardise_columns(coloc)

    if not target_info.empty:
        target_info = standardise_columns(target_info)

    # make protein column consistent before loading into PostgreSQL
    if "protein_id" in mr.columns:
        mr = mr.rename(columns={"protein_id": "protein"})

    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    if not target_info.empty and "protein_id" in target_info.columns:
        target_info = target_info.rename(columns={"protein_id": "protein"})

    # make sure the selected dataset is always recorded
    if "pqtl_dataset" not in mr.columns:
        mr["pqtl_dataset"] = pqtl_dataset

    if "pqtl_dataset" not in coloc.columns:
        coloc["pqtl_dataset"] = pqtl_dataset


    # refresh dashboard tables
    mr.to_sql(mr_table, conn.engine, if_exists="replace", index=False)
    coloc.to_sql(coloc_table, conn.engine, if_exists="replace", index=False)

    phewas_available = not phewas.empty

    if phewas_available:
        phewas.to_sql(phewas_table, conn.engine, if_exists="replace", index=False)

    with st.sidebar.expander("Tracking", expanded=False):
        st.write(f"Loaded {len(mr)} rows into {mr_table}")
        st.write(f"Loaded {len(coloc)} rows into {coloc_table}")


        if phewas_available:
            st.write(f"Loaded {len(phewas)} rows into {phewas_table}")

        if not target_info.empty:
            st.write(f"Loaded {len(target_info)} harmonised top cis-hit rows")

    # load MR + COLOC results
    mr = conn.query(f"SELECT * FROM {mr_table};", ttl=0)
    coloc = conn.query(f"SELECT * FROM {coloc_table};", ttl=0)


    if phewas_available:
        phewas = conn.query(f"SELECT * FROM {phewas_table};", ttl=0)
    else:
        phewas = pd.DataFrame()

    # make PheWAS columns consistent
    if not phewas.empty:
        phewas = phewas.rename(columns={
            "PROTEIN": "protein",
            "protein_id": "protein",
            "PHENO_ID": "pheno_id",
            "OUTCOME_TRAIT": "outcome_trait",
            "PHENOCODE": "phenocode",
            "PHENOSTRING": "phenostring",
            "CATEGORY": "category",
            "SNP": "snp",
            "RSID": "rsid",
            "METHOD": "method",
            "N_INSTRUMENTS": "n_instruments",
            "BETA_MR": "beta_mr",
            "SE_MR": "se_mr",
            "P_MR": "p_mr",
            "P_FDR": "p_fdr",
            "FDR_Q": "fdr_q",
            "FDR_SIGNIFICANT": "fdr_significant",
            "P_BONFERRONI": "p_bonferroni",
            "BONFERRONI_SIGNIFICANT": "bonferroni_significant"
        })

        phewas = standardise_columns(phewas)

        # A1/A2 already come from the original outcome GWAS
        # do not overwrite them with FinnGen ALT/REF
        if "a1" in phewas.columns:
            phewas["a1"] = phewas["a1"].astype(str).str.upper()

        if "a2" in phewas.columns:
            phewas["a2"] = phewas["a2"].astype(str).str.upper()

        for col in [
            "n_instruments",
            "beta_mr",
            "se_mr",
            "p_mr",
            "p_fdr",
            "fdr_q",
            "p_bonferroni"
        ]:
            if col in phewas.columns:
                phewas[col] = pd.to_numeric(phewas[col], errors="coerce")

        if "fdr_significant" in phewas.columns:
            phewas["fdr_significant"] = (
                phewas["fdr_significant"]
                .astype(str)
                .str.lower()
                .isin(["true", "1", "yes"])
            )
        elif "fdr_q" in phewas.columns:
            phewas["fdr_significant"] = phewas["fdr_q"].fillna(np.inf) <= 0.05
        elif "p_fdr" in phewas.columns:
            phewas["fdr_significant"] = phewas["p_fdr"].fillna(np.inf) <= 0.05

        if "bonferroni_significant" in phewas.columns:
            phewas["bonferroni_significant"] = (
                phewas["bonferroni_significant"]
                .astype(str)
                .str.lower()
                .isin(["true", "1", "yes"])
            )
        elif "p_bonferroni" in phewas.columns:
            phewas["bonferroni_significant"] = phewas["p_bonferroni"].fillna(np.inf) <= 0.05

    # MR ammenities
    # standardise numeric MR columns
    mr_numeric_cols = [
        "n_instruments",
        "ivw_beta",
        "ivw_se",
        "ivw_pval",
        "ivw_fdr_q",
        "wald_beta",
        "wald_se",
        "wald_pval",
        "wald_fdr_q",
        "q_pval",
        "egger_intercept_pval"
    ]

    for col in mr_numeric_cols:
        if col not in mr.columns:
            mr[col] = np.nan

        mr[col] = pd.to_numeric(mr[col], errors="coerce")

    coloc_numeric_cols = [
        "pp_h0_abf",
        "pp_h1_abf",
        "pp_h2_abf",
        "pp_h3_abf",
        "pp_h4_abf"
    ]

    for col in coloc_numeric_cols:
        if col in coloc.columns:
            coloc[col] = pd.to_numeric(coloc[col], errors="coerce")

    if "protein_id" in mr.columns:
        mr = mr.rename(columns={"protein_id": "protein"})

    if "protein_id" in coloc.columns:
        coloc = coloc.rename(columns={"protein_id": "protein"})

    if not target_info.empty:
        target_numeric_cols = [
            "frq",
            "gwas_beta",
            "gwas_se",
            "gwas_p",
            "pqtl_beta",
            "pqtl_se",
            "pqtl_p"
        ]

        for col in target_numeric_cols:
            if col in target_info.columns:
                target_info[col] = pd.to_numeric(target_info[col], errors="coerce")

        for col in ["a1", "a2"]:
            if col in target_info.columns:
                target_info[col] = target_info[col].astype(str).str.upper()

        if "protein" in target_info.columns:
            target_info = target_info.drop_duplicates(subset=["protein"])

    required_mr_cols = ["protein", "outcome_trait", "n_instruments"]
    missing_mr_cols = [col for col in required_mr_cols if col not in mr.columns]

    if len(missing_mr_cols) > 0:
        st.error(
            f"cis-MR result file is missing required columns: {missing_mr_cols}. "
            f"File: {mr_file}"
        )
        st.stop()

    required_coloc_cols = ["protein", "outcome_trait"]
    missing_coloc_cols = [col for col in required_coloc_cols if col not in coloc.columns]

    if len(missing_coloc_cols) > 0:
        st.error(
            f"COLOC result file is missing required columns: {missing_coloc_cols}. "
            f"File: {coloc_file}"
        )
        st.stop()

    # if 1 instrument -> use Wald
    # otherwise -> use IVW
    mr["mr_method"] = np.where(mr["n_instruments"] == 1, "Wald ratio", "IVW")
    mr["mr_beta"] = np.where(mr["n_instruments"] == 1, mr["wald_beta"], mr["ivw_beta"])
    mr["mr_se"] = np.where(mr["n_instruments"] == 1, mr["wald_se"], mr["ivw_se"])
    mr["mr_pval"] = np.where(mr["n_instruments"] == 1, mr["wald_pval"], mr["ivw_pval"])
    mr["mr_fdr_q"] = np.where(mr["n_instruments"] == 1, mr["wald_fdr_q"], mr["ivw_fdr_q"])

    # standardise selected pQTL dataset
    selected_pqtl_dataset = (
        str(pqtl_dataset)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )

    # subset database tables to the selected pQTL dataset where possible
    for dataset_col in ["pqtl_dataset", "dataset"]:
        if dataset_col in mr.columns:
            mr_dataset_values = (
                mr[dataset_col]
                .astype(str)
                .str.strip()
                .str.lower()
                .str.replace("-", "_", regex=False)
                .str.replace(" ", "_", regex=False)
            )

            mr = mr[mr_dataset_values == selected_pqtl_dataset].copy()
            break

    for dataset_col in ["pqtl_dataset", "dataset"]:
        if dataset_col in coloc.columns:
            coloc_dataset_values = (
                coloc[dataset_col]
                .astype(str)
                .str.strip()
                .str.lower()
                .str.replace("-", "_", regex=False)
                .str.replace(" ", "_", regex=False)
            )

            coloc = coloc[coloc_dataset_values == selected_pqtl_dataset].copy()
            break

    # available outcomes and default CLI phenotype
    outcomes = sorted(mr["outcome_trait"].dropna().unique())

    if len(outcomes) == 0:
        st.error(f"No cis-MR results were found in {mr_file} for {dataset_name}.")
        st.stop()

    default_outcome = outcomes.index(phenotype) if phenotype in outcomes else 0

    # sidebar filters
    outcome = st.sidebar.selectbox("Outcome", outcomes, index=default_outcome)
    fdr = st.sidebar.slider("MR FDR threshold", 0.0, 1.0, 0.05, 0.01)
    q_pval = st.sidebar.slider("Minimum Cochran Q p-value", 0.0, 1.0, 0.05, 0.01)
    pp4 = st.sidebar.slider("pQTL–GWAS COLOC PP.H4 threshold", 0.0, 1.0, 0.70, 0.01)
    protein = st.sidebar.text_input("Protein search")
    st.title(f"{db_name}: {dataset_name} (N={dataset_n:,}) → {outcome}")

    st.caption(
        f"MR FDR ≤ {fdr:.2f} | Q p ≥ {q_pval:.2f} | "
        f"pQTL–GWAS PP.H4 ≥ {pp4:.2f}"
    )

    # subset everything to selected outcome
    mr_outcome = mr[mr["outcome_trait"] == outcome].copy()
    coloc_outcome = coloc[coloc["outcome_trait"] == outcome].copy()

    if not phewas.empty:
        phewas_outcome = phewas.copy()

        if "outcome_trait" in phewas_outcome.columns:
            phewas_outcome = phewas_outcome[phewas_outcome["outcome_trait"] == outcome].copy()
        elif "pheno_id" in phewas_outcome.columns:
            phewas_outcome = phewas_outcome[phewas_outcome["pheno_id"] == outcome].copy()
    else:
        phewas_outcome = pd.DataFrame()

    # STAGE 1
    # cis-MR supported proteins
    mr_pass = mr_outcome.copy()

    if "mr_fdr_q" in mr_pass.columns:
        mr_pass = mr_pass[mr_pass["mr_fdr_q"].fillna(np.inf) <= fdr]

    # apply Cochran Q only to IVW proteins
    # Wald proteins have no Cochran Q so keep them
    if "q_pval" in mr_pass.columns:
        mr_pass = mr_pass[
            ((mr_pass["mr_method"] == "IVW") & (mr_pass["q_pval"].fillna(-np.inf) >= q_pval))
            |
            (mr_pass["mr_method"] == "Wald ratio")
        ]


    # STAGE 2
    # pQTL - GWAS COLOC
    coloc_pass = coloc_outcome.copy()

    if "pp_h4_abf" not in coloc_pass.columns:
        st.error(
            "The COLOC results do not contain the required PP.H4 column. "
            f"Available columns: {list(coloc_pass.columns)}"
        )
        st.stop()

    coloc_pass["pp_h4_abf"] = pd.to_numeric(
        coloc_pass["pp_h4_abf"],
        errors="coerce"
    )

    coloc_pass = coloc_pass[
        coloc_pass["pp_h4_abf"].fillna(0) >= pp4
    ].copy()

    # proteins which pass both MR + COLOC thresholds
    mr_coloc_pass = mr_pass.merge(
        coloc_pass,
        on="protein",
        how="inner",
        suffixes=("_mr", "_pqtl_coloc")
    )

    # preserve assay-specific protein IDs
    # only remove fully duplicated merged rows
    mr_coloc_pass = mr_coloc_pass.drop_duplicates()

    # add harmonised top cis-hit information
    if not target_info.empty and "protein" in target_info.columns:
        target_cols = [
            "protein",
            "snp",
            "a1",
            "a2",
            "frq",
            "gwas_beta",
            "gwas_se",
            "gwas_p",
            "pqtl_beta",
            "pqtl_se",
            "pqtl_p"
        ]

        target_cols = available_cols(target_info, target_cols)

        mr_coloc_pass = mr_coloc_pass.merge(
            target_info[target_cols],
            on="protein",
            how="left"
        )


    # protein search
    if protein:
        mr_outcome = filter_protein(mr_outcome, protein)
        mr_pass = filter_protein(mr_pass, protein)
        coloc_pass = filter_protein(coloc_pass, protein)
        mr_coloc_pass = filter_protein(mr_coloc_pass, protein)
        phewas_outcome = filter_protein(phewas_outcome, protein)

    # round coloc posterior probs
    for col in coloc_numeric_cols:
        if col in coloc_pass.columns:
            coloc_pass[col] = coloc_pass[col].round(3)

        if col in mr_coloc_pass.columns:
            mr_coloc_pass[col] = mr_coloc_pass[col].round(3)

    # main staged target counts
    n_tested = safe_nunique(mr_outcome, "protein")
    n_mr = safe_nunique(mr_pass, "protein")
    n_mr_coloc = safe_nunique(mr_coloc_pass, "protein")
    n_phewas = safe_nunique(phewas_outcome, "protein")

    tab1, tab2, tab3, tab4 = st.tabs([
        "Overview",
        "cis-MR results",
        "pQTL–GWAS COLOC",
        "Clinical PheWAS"
    ])

    with tab1:
        st.subheader("Target prioritisation")
        metric1, metric2, metric3 = st.columns(3)
        metric1.metric("Proteins tested by cis-MR", n_tested)
        metric2.metric("cis-MR supported", n_mr, f"{retention(n_mr, n_tested):.1f}% of tested", delta_color="off")
        metric3.metric("cis-MR + pQTL COLOC", n_mr_coloc, f"{retention(n_mr_coloc, n_mr):.1f}% retained", delta_color="off")

        funnel_df = pd.DataFrame({
            "stage": [
                "Proteins tested by cis-MR",
                "cis-MR supported",
                "cis-MR + pQTL COLOC"
            ],
            "n_targets": [
                n_tested,
                n_mr,
                n_mr_coloc
            ]
        })

        funnel_fig = px.bar(
            funnel_df,
            x="n_targets",
            y="stage",
            orientation="h",
            text="n_targets",
            title="Progressive target prioritisation",
            labels={"n_targets": "Number of unique proteins", "stage": ""},
            height=420,
            template="plotly_white"
        )

        funnel_fig.update_yaxes(categoryorder="array", categoryarray=funnel_df["stage"][::-1])
        funnel_fig.update_traces(textposition="outside")
        st.plotly_chart(funnel_fig, use_container_width=True)

        st.subheader("Prioritised targets")

        if not mr_coloc_pass.empty:
            prioritised_cols = [
                "protein",
                "snp",
                "a1",
                "a2",
                "frq",
                "gwas_beta",
                "gwas_se",
                "gwas_p",
                "pqtl_beta",
                "pqtl_se",
                "pqtl_p",
                "mr_method",
                "n_instruments",
                "mr_beta",
                "mr_se",
                "mr_pval",
                "mr_fdr_q",
                "q_pval",
                "egger_intercept_pval",
                "pp_h4_abf"
            ]

            prioritised_cols = available_cols(mr_coloc_pass, prioritised_cols)

            if "pp_h4_abf" in mr_coloc_pass.columns:
                mr_coloc_pass = mr_coloc_pass.sort_values(
                    ["pp_h4_abf", "mr_fdr_q"],
                    ascending=[False, True],
                    na_position="last"
                )

            overview_table = mr_coloc_pass[prioritised_cols].copy()

            overview_column_names = {
                "protein": "Protein",
                "snp": "Top SNP",
                "a1": "Risk allele",
                "a2": "Other allele",
                "frq": "Risk allele frequency",
                "gwas_beta": "GWAS beta",
                "gwas_se": "GWAS SE",
                "gwas_p": "GWAS p-value",
                "pqtl_beta": "pQTL beta",
                "pqtl_se": "pQTL SE",
                "pqtl_p": "pQTL p-value",
                "mr_method": "MR method",
                "n_instruments": "N instruments",
                "mr_beta": "MR beta",
                "mr_se": "MR SE",
                "mr_pval": "MR p-value",
                "mr_fdr_q": "MR FDR",
                "q_pval": "Cochran Q p-value",
                "egger_intercept_pval": "Egger intercept p-value",
                "pp_h4_abf": "COLOC PP.H4"
            }

            overview_table = overview_table.rename(columns=overview_column_names)

            st.success(
                f"{n_mr_coloc} unique target(s) passed the selected cis-MR and pairwise COLOC thresholds."
            )

            st.caption(
                "All SNP effects are harmonised to the outcome GWAS risk allele. "
                "A positive pQTL beta means the risk allele increases protein abundance, "
                "whereas a negative pQTL beta means the risk allele decreases protein abundance."
            )

            st.dataframe(
                overview_table,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Risk allele frequency": st.column_config.NumberColumn(format="%.3f"),
                    "GWAS beta": st.column_config.NumberColumn(format="%.4f"),
                    "GWAS SE": st.column_config.NumberColumn(format="%.4f"),
                    "GWAS p-value": st.column_config.NumberColumn(format="%.3e"),
                    "pQTL beta": st.column_config.NumberColumn(format="%.4f"),
                    "pQTL SE": st.column_config.NumberColumn(format="%.4f"),
                    "pQTL p-value": st.column_config.NumberColumn(format="%.3e"),
                    "MR beta": st.column_config.NumberColumn(format="%.4f"),
                    "MR SE": st.column_config.NumberColumn(format="%.4f"),
                    "MR p-value": st.column_config.NumberColumn(format="%.3e"),
                    "MR FDR": st.column_config.NumberColumn(format="%.3e"),
                    "Cochran Q p-value": st.column_config.NumberColumn(format="%.3e"),
                    "Egger intercept p-value": st.column_config.NumberColumn(format="%.3e"),
                    "COLOC PP.H4": st.column_config.NumberColumn(format="%.3f")
                }
            )

            st.download_button(
                label="Download prioritised targets",
                data=overview_table.to_csv(index=False, sep="\t"),
                file_name=f"{pqtl_dataset}_{outcome}_prioritised_target_overview.tsv",
                mime="text/tab-separated-values",
                key="download_prioritised_targets_overview"
            )

        else:
            st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

        if phewas_available:
            st.caption(
                f"PheWAS safety results are available for {n_phewas} unique target(s)."
            )

    with tab2:
        show_all_mr = st.checkbox("Show all tested cis-MR proteins", value=False)
        mr_display = mr_outcome if show_all_mr else mr_pass
        n_ivw = (mr_display["mr_method"] == "IVW").sum()
        n_wald = (mr_display["mr_method"] == "Wald ratio").sum()

        col1, col2, col3 = st.columns(3)
        col1.metric("MR proteins shown", mr_display["protein"].nunique())
        col2.metric("IVW proteins", int(n_ivw))
        col3.metric("Wald proteins", int(n_wald))

        display_cols = [
            "protein",
            "outcome_trait",
            "n_instruments",
            "mr_method",
            "mr_beta",
            "mr_se",
            "mr_pval",
            "mr_fdr_q",
            "q_pval",
            "egger_intercept_pval",
            "ivw_beta",
            "ivw_se",
            "ivw_pval",
            "ivw_fdr_q",
            "wald_beta",
            "wald_se",
            "wald_pval",
            "wald_fdr_q"
        ]

        display_cols = available_cols(mr_display, display_cols)
        remaining_cols = [col for col in mr_display.columns if col not in display_cols]

        st.dataframe(
            mr_display[display_cols + remaining_cols],
            use_container_width=True,
            hide_index=True
        )

        # primary MR volcano plot
        plot_df = mr_display[
            mr_display["mr_pval"].notna() &
            mr_display["mr_beta"].notna() &
            (mr_display["mr_pval"] > 0)
        ].copy()

        if not plot_df.empty:
            plot_df["minus_log10_mr_pval"] = -np.log10(plot_df["mr_pval"])
            plot_df["significant"] = plot_df["mr_fdr_q"] < 0.05

            fig = px.scatter(
                plot_df,
                x="mr_beta",
                y="minus_log10_mr_pval",
                hover_name="protein",
                color="significant",
                symbol="mr_method",
                hover_data={
                    "mr_method": True,
                    "n_instruments": True,
                    "mr_beta": ":.4f",
                    "mr_se": ":.4f",
                    "mr_pval": ":.3e",
                    "mr_fdr_q": ":.3e",
                    "minus_log10_mr_pval": False
                },
                labels={
                    "mr_beta": "Primary MR beta",
                    "minus_log10_mr_pval": "-log10(primary MR p-value)",
                    "mr_method": "MR method",
                    "significant": "FDR < 0.05"
                },
                title="Primary cis-MR volcano plot",
                height=600,
                template="plotly_white"
            )

            fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey")
            fig.add_vline(x=0, line_dash="dash", line_color="grey")
            st.plotly_chart(fig, use_container_width=True)

        else:
            st.info("No MR results remain after applying the selected filters.")

    with tab3:
        st.subheader("cis-MR + pQTL–GWAS COLOC targets")
        st.caption("Targets shown here pass both the selected cis-MR and pairwise pQTL–GWAS COLOC thresholds.")

        if not mr_coloc_pass.empty:
            col1, col2, col3 = st.columns(3)
            col1.metric("Prioritised proteins", mr_coloc_pass["protein"].nunique())
            col2.metric("Median PP.H4", f"{mr_coloc_pass['pp_h4_abf'].median():.3f}" if "pp_h4_abf" in mr_coloc_pass.columns else "NA")
            col3.metric("Median MR FDR", f"{mr_coloc_pass['mr_fdr_q'].median():.3e}" if "mr_fdr_q" in mr_coloc_pass.columns else "NA")

            prioritised_cols = [
                "protein",
                "mr_method",
                "n_instruments",
                "mr_beta",
                "mr_se",
                "mr_pval",
                "mr_fdr_q",
                "q_pval",
                "egger_intercept_pval",
                "top_snp",
                "pp_h0_abf",
                "pp_h1_abf",
                "pp_h2_abf",
                "pp_h3_abf",
                "pp_h4_abf"
            ]

            prioritised_cols = available_cols(mr_coloc_pass, prioritised_cols)
            remaining_cols = [col for col in mr_coloc_pass.columns if col not in prioritised_cols]

            if "pp_h4_abf" in mr_coloc_pass.columns:
                mr_coloc_pass = mr_coloc_pass.sort_values(
                    ["pp_h4_abf", "mr_fdr_q"],
                    ascending=[False, True],
                    na_position="last"
                )

            st.dataframe(
                mr_coloc_pass[prioritised_cols + remaining_cols],
                use_container_width=True,
                hide_index=True
            )

            st.download_button(
                label="Download prioritised targets",
                data=mr_coloc_pass.to_csv(index=False, sep="\t"),
                file_name=f"{outcome}_prioritised_targets.tsv",
                mime="text/tab-separated-values",
                key="download_prioritised_targets_coloc"
            )

        else:
            st.info("No proteins currently pass both the selected cis-MR and pQTL COLOC thresholds.")

    with tab4:
        st.subheader("FinnGen PheWAS safety and repurposing profile")

        if phewas_outcome.empty:
            st.info("No local PheWAS safety results are available for this outcome.")

        elif "protein" not in phewas_outcome.columns:
            st.error("The PheWAS result file does not contain a protein column.")

        else:
            phewas_targets = sorted(phewas_outcome["protein"].dropna().astype(str).unique())

            if len(phewas_targets) == 0:
                st.info("No proteins were found in the PheWAS safety table.")

            else:
                default_phewas_target = 0
                prioritised_target_names = mr_coloc_pass["protein"].dropna().astype(str).unique().tolist()

                for target in prioritised_target_names:
                    if target in phewas_targets:
                        default_phewas_target = phewas_targets.index(target)
                        break

                selected_phewas_target = st.selectbox(
                    "Target",
                    phewas_targets,
                    index=default_phewas_target,
                    key="selected_phewas_target"
                )

                target_phewas = phewas_outcome[
                    phewas_outcome["protein"].astype(str) == selected_phewas_target
                ].copy()

                p_col = None
                beta_col = None
                bonferroni_col = None

                for col in ["p_mr"]:
                    if col in target_phewas.columns:
                        p_col = col
                        break

                for col in ["beta_mr"]:
                    if col in target_phewas.columns:
                        beta_col = col
                        break

                for col in ["p_bonferroni"]:
                    if col in target_phewas.columns:
                        bonferroni_col = col
                        break

                if p_col is None or beta_col is None:
                    st.error(
                        "The PheWAS result file needs the MR effect column "
                        "(beta_mr) and the MR p-value column "
                        "(p_mr)."
                    )

                else:
                    target_phewas = target_phewas[
                        target_phewas[p_col].notna() &
                        target_phewas[beta_col].notna() &
                        (target_phewas[p_col] > 0)
                    ].copy()

                    if target_phewas.empty:
                        st.info(f"No valid PheWAS associations were found for {selected_phewas_target}.")

                    else:
                        target_phewas["minus_log10_p"] = -np.log10(target_phewas[p_col])

                        if bonferroni_col is not None:
                            target_phewas["bonferroni_significant"] = target_phewas[bonferroni_col].fillna(np.inf) <= 0.05
                        elif "bonferroni_significant" not in target_phewas.columns:
                            target_phewas["bonferroni_significant"] = False

                        phenotype_col = "phenostring" if "phenostring" in target_phewas.columns else "phenocode"
                        category_col = "category" if "category" in target_phewas.columns else None

                        n_phenotypes = target_phewas[phenotype_col].nunique()
                        n_nominal = int((target_phewas[p_col] < 0.05).sum())
                        n_bonferroni = int(target_phewas["bonferroni_significant"].sum())

                        metric1, metric2, metric3 = st.columns(3)
                        metric1.metric("FinnGen phenotypes tested", int(n_phenotypes))
                        metric2.metric("Nominal associations", n_nominal)
                        metric3.metric("Bonferroni-significant associations", n_bonferroni)
                        
                        st.caption(
                            "PheWAS MR estimates show the effect of genetically predicted protein abundance "
                            "on each FinnGen phenotype (ICD coded). Wald ratio is used for targets with "
                            "one available cis-MR instrument and IVW is used for targets with > 1."
                        )

                        plot_kwargs = {
                            "data_frame": target_phewas,
                            "x": beta_col,
                            "y": "minus_log10_p",
                            "hover_name": phenotype_col,
                            "symbol": "bonferroni_significant",
                            "hover_data": {
                                beta_col: ":.4f",
                                p_col: ":.3e",
                                "minus_log10_p": False,
                                "bonferroni_significant": True
                            },
                            "labels": {
                                beta_col: "PheWAS MR beta",
                                "minus_log10_p": "-log10(PheWAS p-value)",
                                "bonferroni_significant": "Bonferroni significant"
                            },
                            "title": f"FinnGen PheWAS profile: {selected_phewas_target}",
                            "height": 600
                        }

                        if "phenocode" in target_phewas.columns:
                            plot_kwargs["hover_data"]["phenocode"] = True

                        if category_col is not None:
                            plot_kwargs["color"] = category_col
                            plot_kwargs["labels"][category_col] = "FinnGen category"

                        phewas_fig = px.scatter(**plot_kwargs)
                        phewas_fig.add_hline(y=-np.log10(0.05 / 2511), line_dash="dash", line_color="grey")
                        phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
                        st.plotly_chart(phewas_fig, use_container_width=True)

                        st.subheader("Bonferroni-significant PheWAS associations")

                        top_phewas = target_phewas[target_phewas["bonferroni_significant"]].copy()

                        if bonferroni_col is not None:
                            top_phewas = top_phewas.sort_values(bonferroni_col, ascending=True)
                        else:
                            top_phewas = top_phewas.sort_values(p_col, ascending=True)

                        top_phewas = top_phewas.sort_values(beta_col, ascending=True)

                        if top_phewas.empty:
                            st.info(
                                f"No FinnGen phenotype associations survive Bonferroni correction across "
                                f"2,511 ICD endpoints for {selected_phewas_target}."
                            )

                        else:
                            top_plot_kwargs = {
                                "data_frame": top_phewas,
                                "x": beta_col,
                                "y": phenotype_col,
                                "hover_data": {
                                    beta_col: ":.4f",
                                    p_col: ":.3e",
                                    "minus_log10_p": ":.3f"
                                },
                                "labels": {
                                    beta_col: "PheWAS MR beta",
                                    phenotype_col: ""
                                },
                                "title": "Bonferroni-significant PheWAS associations",
                                "height": max(450, 45 * len(top_phewas))
                            }

                            if "phenocode" in top_phewas.columns:
                                top_plot_kwargs["hover_data"]["phenocode"] = True

                            if category_col is not None:
                                top_plot_kwargs["color"] = category_col
                                top_plot_kwargs["labels"][category_col] = "FinnGen category"

                            top_phewas_fig = px.scatter(**top_plot_kwargs)
                            top_phewas_fig.add_vline(x=0, line_dash="dash", line_color="grey")
                            st.plotly_chart(top_phewas_fig, use_container_width=True)

                        phewas_cols = [
                            "protein",
                            "method",
                            "n_instruments",
                            "rsid",
                            "A1",
                            "A2",
                            "phenocode",
                            "phenostring",
                            "category",
                            "beta_mr",
                            "se_mr",
                            "p_mr",
                            "p_bonferroni",
                            "bonferroni_significant"
                        ]

                        phewas_cols = [
                            col for col in phewas_cols
                            if col is not None and col in target_phewas.columns
                        ]

                        significant_phewas = target_phewas[
                            target_phewas["bonferroni_significant"]
                        ].copy()

                        if bonferroni_col is not None:
                            significant_phewas = significant_phewas.sort_values(bonferroni_col, ascending=True)
                        else:
                            significant_phewas = significant_phewas.sort_values(p_col, ascending=True)

                        if significant_phewas.empty:
                            st.success("No FinnGen phenotype associations survive Bonferroni correction across 2,511 ICD endpoints for this target.")
                        else:
                            st.dataframe(
                                significant_phewas[phewas_cols],
                                use_container_width=True,
                                hide_index=True
                            )

                        with st.expander("View all PheWAS associations"):
                            remaining_cols = [col for col in target_phewas.columns if col not in phewas_cols]
                            st.dataframe(
                                target_phewas[phewas_cols + remaining_cols].sort_values(p_col, ascending=True),
                                use_container_width=True,
                                hide_index=True
                            )

                        st.download_button(
                            label=f"Download {selected_phewas_target} PheWAS results",
                            data=target_phewas.to_csv(index=False, sep="\t"),
                            file_name=f"{selected_phewas_target}_{outcome}_FinnGen_PheWAS.tsv",
                            mime="text/tab-separated-values",
                            key=f"download_phewas_{pqtl_dataset}_{outcome}_{selected_phewas_target}"
                        )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db_name", required=True, type=str)
    p.add_argument("--port_number", required=True, type=str)
    p.add_argument("--phenotype", required=True, type=str)
    p.add_argument("--pqtl_dataset", required=True, type=str)
    args = p.parse_args()
    create_streamlit_ammenities(args.db_name, args.port_number)
    dashboard(
        db_name=args.db_name,
        port_number=args.port_number,
        phenotype=args.phenotype,
        pqtl_dataset=args.pqtl_dataset
    )


if __name__ == "__main__":
    main()