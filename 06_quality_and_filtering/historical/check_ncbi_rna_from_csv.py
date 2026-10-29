#!/usr/bin/env python3
import csv
import re
import time
import requests
import xml.etree.ElementTree as ET
from pathlib import Path

INPUT = Path.home() / "Desktop" / "rna_check.csv"
OUTPUT = Path.home() / "Desktop" / "rna_sra_summary.csv"

HEADERS = {"User-Agent": "Mozilla/5.0"}
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"

def fetch_text(url, params=None, timeout=120):
    r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    return r.text

def esearch(db, term, retmax=20):
    txt = fetch_text(BASE + "esearch.fcgi", {
        "db": db,
        "term": term,
        "retmode": "xml",
        "retmax": str(retmax)
    })
    root = ET.fromstring(txt)
    return [x.text for x in root.findall(".//Id") if x.text]

def esummary_assembly(uid):
    txt = fetch_text(BASE + "esummary.fcgi", {
        "db": "assembly",
        "id": uid,
        "retmode": "xml"
    })
    root = ET.fromstring(txt)
    info = {
        "assembly": "",
        "species": "",
        "biosample": "",
        "bioproject": ""
    }
    for doc in root.findall(".//DocumentSummary"):
        info["assembly"] = doc.findtext("AssemblyAccession") or ""
        info["species"] = doc.findtext("Organism") or ""
        info["biosample"] = doc.findtext("BioSampleAccn") or ""
        info["bioproject"] = doc.findtext("BioProjectAccn") or ""
    return info

def elink_ids(dbfrom, db, uid):
    txt = fetch_text(BASE + "elink.fcgi", {
        "dbfrom": dbfrom,
        "db": db,
        "id": uid,
        "retmode": "xml"
    })
    root = ET.fromstring(txt)
    return sorted(set(x.text for x in root.findall(".//LinkSetDb/Link/Id") if x.text))

def get_runinfo_by_acc(acc):
    try:
        txt = fetch_text(RUNINFO_URL, {"acc": acc}, timeout=120)
        return txt
    except Exception:
        return ""

def parse_runinfo(text):
    if not text.strip():
        return []
    lines = text.splitlines()
    if len(lines) < 2:
        return []
    return list(csv.DictReader(lines))

def clean_acc(x):
    x = (x or "").strip()
    if re.match(r"^GC[AF]_\d+\.\d+$", x):
        return x
    return ""

def uniq_join(vals):
    vals = [v for v in vals if v]
    vals = sorted(set(vals))
    return ";".join(vals)

def is_rna_row(row):
    joined = " ".join([
        row.get("LibraryStrategy", "") or "",
        row.get("LibrarySource", "") or "",
        row.get("LibrarySelection", "") or "",
        row.get("ScientificName", "") or "",
        row.get("SampleName", "") or "",
        row.get("StudyTitle", "") or "",
        row.get("Experiment", "") or ""
    ]).lower()

    keywords = ["rna-seq", "rna seq", "rnaseq", "transcript", "mrna", "total rna"]
    return any(k in joined for k in keywords)

def is_usable_row(row):
    strategy = (row.get("LibraryStrategy", "") or "").strip().upper()
    source = (row.get("LibrarySource", "") or "").strip().lower()
    platform = (row.get("Platform", "") or "").strip().upper()

    joined = " ".join([
        row.get("LibraryStrategy", "") or "",
        row.get("LibrarySource", "") or "",
        row.get("LibrarySelection", "") or "",
        row.get("StudyTitle", "") or "",
        row.get("Experiment", "") or "",
        row.get("SampleName", "") or ""
    ]).lower()

    bad_keywords = [
        "small rna", "mirna", "srna", "amplicon", "metatranscript",
        "single cell", "chip-seq", "atac", "bisulfite"
    ]
    if any(k in joined for k in bad_keywords):
        return False

    if strategy == "RNA-SEQ" and source in ["transcriptomic", "transcriptome", ""]:
        return True

    if "transcript" in source and platform in ["ILLUMINA", ""]:
        return True

    return False

records = []
with open(INPUT, "r", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f)
    reader.fieldnames = [x.strip() if x else x for x in reader.fieldnames]

    for row in reader:
        row = {(k.strip() if k else k): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
        sample = (row.get("sample") or "").strip()
        acc = clean_acc(row.get("GCA_ID") or "")
        if acc:
            records.append({"sample": sample, "assembly": acc})

out_rows = []

for n, rec in enumerate(records, 1):
    acc = rec["assembly"]
    sample = rec["sample"]

    result = {
        "assembly": acc,
        "species": sample,
        "biosample": "",
        "bioproject": "",
        "linked_sra_total": 0,
        "rna_sra_records": 0,
        "usable_runs": "",
        "library_strategy": "",
        "platform": "",
        "note": ""
    }

    try:
        assembly_ids = esearch("assembly", acc, retmax=5)
        time.sleep(0.34)

        if not assembly_ids:
            result["note"] = "assembly_not_found_in_ncbi"
            out_rows.append(result)
            continue

        assembly_uid = assembly_ids[0]
        meta = esummary_assembly(assembly_uid)
        time.sleep(0.34)

        result["species"] = meta["species"] or sample
        result["biosample"] = meta["biosample"]
        result["bioproject"] = meta["bioproject"]

        sra_ids = elink_ids("assembly", "sra", assembly_uid)
        time.sleep(0.34)

        if not sra_ids and meta["biosample"]:
            biosample_ids = esearch("biosample", meta["biosample"], retmax=5)
            time.sleep(0.34)
            if biosample_ids:
                sra_ids = elink_ids("biosample", "sra", biosample_ids[0])
                time.sleep(0.34)

        result["linked_sra_total"] = len(sra_ids)

        run_rows = []
        if meta["bioproject"]:
            run_rows = parse_runinfo(get_runinfo_by_acc(meta["bioproject"]))

        if not run_rows and meta["biosample"]:
            run_rows = parse_runinfo(get_runinfo_by_acc(meta["biosample"]))

        rna_rows = [r for r in run_rows if is_rna_row(r)]
        usable_rows = [r for r in rna_rows if is_usable_row(r)]

        result["rna_sra_records"] = len(rna_rows)
        result["usable_runs"] = uniq_join([r.get("Run", "") for r in usable_rows])
        result["library_strategy"] = uniq_join([r.get("LibraryStrategy", "") for r in usable_rows])
        result["platform"] = uniq_join([r.get("Platform", "") for r in usable_rows])

        if result["rna_sra_records"] == 0:
            result["note"] = "no_rna_found_via_bioproject_or_biosample"
        elif len(usable_rows) == 0:
            result["note"] = "rna_found_but_not_recommended_for_braker3"
        else:
            result["note"] = "usable_rna_found"

    except Exception as e:
        result["note"] = f"error:{str(e)[:160].replace(',', ';')}"

    out_rows.append(result)
    print(f"[{n}/{len(records)}] {acc} done")

with open(OUTPUT, "w", encoding="utf-8-sig", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "assembly", "species", "biosample", "bioproject",
        "linked_sra_total", "rna_sra_records", "usable_runs",
        "library_strategy", "platform", "note"
    ])
    writer.writeheader()
    writer.writerows(out_rows)

print(f"Done. Output: {OUTPUT}")
print(f"Processed NCBI accessions: {len(out_rows)}")
