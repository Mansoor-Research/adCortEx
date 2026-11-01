"""Freeze sample and probe metadata before any external association analysis."""

import csv
import gzip
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/external/GSE33000/GSE33000_series_matrix.txt.gz"
ANNOTATION = ROOT / "data/external/GSE33000/GPL4372.annot.gz"
OUT = ROOT / "results/external_validation"
OUT.mkdir(parents=True, exist_ok=True)

fields = {}
with gzip.open(SOURCE, "rt", encoding="utf-8", errors="replace", newline="") as handle:
    for row in csv.reader(handle, delimiter="\t"):
        if row and row[0] == "!series_matrix_table_begin":
            break
        if not row or not row[0].startswith("!Sample_"):
            continue
        fields.setdefault(row[0], []).append(row[1:])

accessions = fields["!Sample_geo_accession"][0]
n = len(accessions)
characteristics = fields["!Sample_characteristics_ch2"]
if any(len(values) != n for values in characteristics):
    raise RuntimeError("Unequal sample metadata lengths")

sample_rows = []
for i, accession in enumerate(accessions):
    attributes = {}
    for line in characteristics:
        key, separator, value = line[i].partition(": ")
        if separator:
            attributes[key.lower()] = value
    age_text = attributes.get("age", "")
    age = age_text.removesuffix(" yrs") if age_text.endswith(" yrs") else ""
    sample_rows.append({
        "gsm": accession,
        "disease_status": attributes.get("disease status", ""),
        "sex": attributes.get("gender", ""),
        "age_years": age,
    })

with (OUT / "sample_metadata.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(sample_rows[0]))
    writer.writeheader()
    writer.writerows(sample_rows)

probe_rows = []
with gzip.open(ANNOTATION, "rt", encoding="utf-8", errors="replace", newline="") as handle:
    for line in handle:
        if line.rstrip() == "!platform_table_begin":
            break
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["ID"] == "!platform_table_end":
            break
        raw_entrez = row.get("Gene ID", "").strip()
        raw_symbol = row.get("Gene symbol", "").strip()
        valid = raw_entrez.isdecimal() and "///" not in raw_symbol and bool(raw_symbol)
        probe_rows.append({
            "probe_id": row["ID"],
            "entrez_id": raw_entrez if valid else "",
            "symbol": raw_symbol if valid else "",
            "mapping_status": "unique" if valid else "excluded_ambiguous_or_missing",
        })

with (OUT / "probe_mapping.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(probe_rows[0]))
    writer.writeheader()
    writer.writerows(probe_rows)

qc = {
    "n_samples": n,
    "disease_counts": dict(Counter(row["disease_status"] for row in sample_rows)),
    "sex_counts": dict(Counter(row["sex"] for row in sample_rows)),
    "age_complete": sum(bool(row["age_years"]) for row in sample_rows),
    "n_probes_annotated": len(probe_rows),
    "n_probes_uniquely_mapped": sum(row["mapping_status"] == "unique" for row in probe_rows),
    "n_unique_entrez": len({row["entrez_id"] for row in probe_rows if row["entrez_id"]}),
}
(OUT / "metadata_qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
print(json.dumps(qc, indent=2))
