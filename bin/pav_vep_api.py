import requests

r = requests.post(
    "https://rest.ensembl.org/vep/human/id",
    headers={"Content-Type": "application/json"},
    json={"ids": ["rs1859788", "rs5848"], "mane": 1, "hgvs": 1, "protein": 1},
)

for v in r.json():
    for tc in v.get("transcript_consequences", []):
        if tc.get("mane_select"):
            print(v["input"], tc["gene_symbol"],
                  tc["consequence_terms"], tc.get("hgvsp"))
