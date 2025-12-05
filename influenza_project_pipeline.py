
# influenza_project_pipeline.py
# Complete pipeline for HA protein sequence processing, metadata extraction,
# year-splitting, and ML feature engineering for influenza evolution analysis.

from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import os

# -------------------------
# LOAD FASTA SEQUENCES
# -------------------------

FASTA_PATH = "data/FASTA.fa"  # replace with your path
records = list(SeqIO.parse(FASTA_PATH, "fasta"))

print(f"Loaded {len(records)} sequences.")

# -------------------------
# METADATA EXTRACTION
# -------------------------

meta_rows = []

for r in records:
    desc_parts = r.description.split()
    accession = desc_parts[0]
    strain = desc_parts[1] if len(desc_parts) > 1 else None

    # Extract year from second-last token (format: 2002// or 2011/01/26)
    year_token = desc_parts[-2]
    m = re.search(r"(\d{4})", year_token)
    year = int(m.group(1)) if m else None

    protein = desc_parts[-1]
    length = len(r.seq)

    meta_rows.append([accession, strain, year, protein, length])

df_meta = pd.DataFrame(meta_rows, columns=["Accession", "Strain", "Year", "Protein", "Length"])

# Save metadata
df_meta.to_csv("metadata.csv", index=False)
print("Metadata saved to metadata.csv")

# -------------------------
# SPLIT SEQUENCES BY YEAR (2010–PRESENT)
# -------------------------

df_recent = df_meta[df_meta["Year"] >= 2010]

output_dir = "HA_by_year"
os.makedirs(output_dir, exist_ok=True)

record_dict = {rec.id: rec for rec in records}

for year in sorted(df_recent["Year"].unique()):
    year_records = df_recent[df_recent["Year"] == year]["Accession"]
    out_path = os.path.join(output_dir, f"HA_{year}.fasta")

    with open(out_path, "w") as f:
        for acc in year_records:
            SeqIO.write(record_dict[acc], f, "fasta")

    print(f"Saved {len(year_records)} sequences for {year} → {out_path}")

# -------------------------
# FEATURE ENGINEERING (AA composition)
# -------------------------

AAS = "ACDEFGHIKLMNPQRSTVWY"

def aa_composition(seq):
    total = len(seq)
    return {aa: seq.count(aa) / total for aa in AAS}

feature_rows = []

for r in records:
    seq = str(r.seq)
    comp = aa_composition(seq)

    desc_parts = r.description.split()
    year_token = desc_parts[-2]
    m = re.search(r"(\d{4})", year_token)
    year = int(m.group(1)) if m else None

    row = {"Accession": r.id, "Year": year, "Length": len(seq)}
    row.update(comp)

    feature_rows.append(row)

df_features = pd.DataFrame(feature_rows)
df_features.to_csv("HA_features.csv", index=False)

print("Feature dataset saved to HA_features.csv")

print("Pipeline completed successfully.")
