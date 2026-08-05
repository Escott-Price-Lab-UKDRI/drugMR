# drugMR: A Multi-Fluid Multi-Omics Drug Discovery Pipeline

End-to-end local + cloud hybrid pipeline leveraging plasma, CSF and brain protein QTLs (pQTLs; >10,000 proteins) from UKBB-PPP, deCODE, WU, and Wingo et al. (Olink + SomaScan + MassSpec) for any phenotype of interest (demonstrated in Alzheimer's disease), with or without mediating biomarkers. From raw GWAS and pQTL data to single-cell target prioritisation, mediation analysis and drug safety assessment. Outputs are compiled into a production-ready dashboard.

---

## Pipeline overview

Each stage's output gates the next (a stage is skipped and reused if its output already exists, unless `overwrite: true`):

1. **GWAS QC** (`bin/qc_gwas.py`) - harmonises and QCs your outcome GWAS.
2. **Mediator QC** (optional, `mediators: true`) - QCs mediating biomarker GWAS.
3. **cis-region preparation** (`bin/prep_cis_regions.py`) - matches pQTL cis-regions to the outcome GWAS.
4. **cis-MR** (`bin/cis_mr.R`) — Wald ratio (1 instrument) / IVW (>1 instrument) MR per protein.
5. **NetworkMR** (optional, with mediators) - mediation analysis through the specified biomarkers.
6. **Pairwise COLOC** (`bin/coloc_targets.py`) - pQTL–GWAS colocalisation (PP.H4) per protein.
7. **Top cis-hit compilation** (`bin/compile_cis_hit_info.py`) - harmonises the top cis-SNP per protein, aligning alleles to the outcome risk allele.
8. **SMR** (`bin/sort_smr.py`) - for targets which survive cis-MR + COLOC, tests eQTL–GWAS colocalisation via SMR + HEIDI, both **bulk** (pre-computed eQTLGen / MetaBrain / GTEx_v10, ingested as-is) and **single-cell** (SingleBrain, computed fresh per cell type). FDR-corrected per dataset, alleles aligned to the outcome risk allele, and for single-cell targets the eQTL beta is re-sourced from the original per-cell-type eQTL file rather than the raw SMR output.
9. **PheWAS safety screening** (`bin/phewas_cis_pqtls.py`, `bin/ukb_phewas.py`) - FinnGen and UKB Biobank phenome-wide MR for surviving targets.
10. **Results** (`dm.results()`) - loads cis-MR/COLOC results into PostgreSQL and launches the Streamlit dashboard.

---

## Clone the repo!

```bash
git clone https://github.com/guillermocomesanacimadevila/drugMR.git
cd drugMR/
```

## Configure `assets/config.yaml`

Every run is driven by `assets/config.yaml`. Key fields:

| Field | Purpose |
| --- | --- |
| `pheno_id`, `sumstats`, `n_cases`, `n_controls` | Outcome GWAS identity and sample size |
| `snp_col` / `a1_col` / `a2_col` / `beta_col` / `se_col` / `p_col` / `pos_col` / `chr_col` / `af_col` | Column names in your outcome GWAS |
| `pqtl_dataset`, `pqtl_dir` | Which pQTL dataset to run (`ukb_ppp`, `decode`, `wu_csf`, `wingo_brain`) |
| `ref_bfile` | Reference panel (1000 Genomes) for cis-MR / SMR |
| `mediators`, `mediator_manifest` | Enable mediation analysis through a manifest of biomarkers |
| `run_smr` | Master on/off switch for the SMR step |
| `bulk_eqtl_datasets` | Pre-computed bulk eQTL datasets to ingest (e.g. `[eQTLGen, MetaBrain, GTEx_v10]`); `[]` skips bulk SMR |
| `sc_eqtl_dataset` | Single-cell eQTL dataset to run SMR against (e.g. `SingleBrain`); empty skips single-cell SMR |
| `maf`, `remove_mhc`, `remove_apoe` | QC filters applied to GWAS/pQTLs |
| `overwrite` | Force every stage to rerun instead of reusing existing outputs |

## Running the pipeline

From a notebook or script:

```python
import drugmr as dm

# run locally via Docker
dm.local(config="assets/config.yaml")

# OR run on the Falcon HPC cluster via SLURM/apptainer
dm.hpc(config="assets/config.yaml")

# load cis-MR/COLOC results into PostgreSQL and launch the Streamlit dashboard
dm.results()
```

See `notebooks/00_drugmr.ipynb` for a worked example.

## Synapse configuration

Create the Synapse config file:

```bash
nano ~/.synapseConfig
```

Populate it with your Synapse info!

```bash
[default]
username = your_email@example.com
authtoken = YOUR_PERSONAL_ACCESS_TOKEN

[cache]
location = ~/.synapseCache
```

## Streamlit configuration

Create the Streamlit secrets file:

```bash
nano .streamlit/secrets.toml
```

Populate `.streamlit/secrets.toml` as follows:

```toml
[connections.postgresql]
dialect = "postgresql"
host = "localhost"
port = "5432"
database = "xxx"
username = "xxx"
password = "xxx"
```
## Configure passwordless SSH access to Falcon

Generate an SSH key (if you do not already have one):

```bash
ssh-keygen -t ed25519 -C "drugMR"
```

Display your public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the public key to Falcon:

```bash
ssh-copy-id c.<username>@falconlogin.cf.ac.uk
```

Test the connection:

```bash
ssh c.<username>@falconlogin.cf.ac.uk
```

## Dashboard

`dm.results()` launches a Streamlit dashboard (`dashboard/mr_app.py`) with a pQTL dataset selector and sidebar filters (outcome, FDR/Q/PP.H4 thresholds, protein search) shared across:

- **Overview** - target prioritisation funnel (cis-MR → cis-MR + COLOC) and the prioritised target table.
- **1. cis-MR** - full cis-MR association table and a volcano plot.
- **2. pQTL–GWAS COLOC** - targets passing both cis-MR and pairwise COLOC thresholds.
- **3. FinnGen PheWAS** / **4. UKB PheWAS** - per-target phenome-wide MR safety profile, with a Manhattan-style scatter and Bonferroni-significant associations.
- **5. SMR (bulk/sc eQTL)** - filterable SMR + HEIDI results (by data type / cell type), plus a target × cell-type/tissue support heatmap.
- **6. Final Targets** - the curated, filter-free deliverable table: targets passing cis-MR + COLOC + SMR + HEIDI, one row per target × cell-type/tissue, with GWAS/pQTL/eQTL/SMR betas all aligned to the outcome risk allele, the top SMR SNP (with chromosome/position), and SMR/HEIDI p-values.

## Docker 

**This is simply a dev note, do NOT worry about it**

Pull the latest `drugMR` image from GHCR:

```bash
docker pull ghcr.io/guillermocomesanacimadevila/drugmr:latest
```

## Authors

**Guillermo Comesaña Cimadevila**<sup>1,2</sup>, **Marie-Joe Dib**<sup>3</sup>, **Dervis Salih**<sup>4</sup>, **Nicholas J. Bray**<sup>2</sup>, **Emily Simmonds**<sup>1</sup>, **Valentina Escott-Price**<sup>1,2</sup>

- <sup>1</sup> UK Dementia Research Institute at Cardiff University, Cardiff, UK.
- <sup>2</sup> MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University, Cardiff, UK.
- <sup>3</sup> Nascent Studio Ltd, London, UK.
- <sup>4</sup> UK Dementia Research Institute at University College London, London, UK.
