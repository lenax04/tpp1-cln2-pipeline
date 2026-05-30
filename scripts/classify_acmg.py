#!/usr/bin/env python3
"""
scripts/classify_acmg.py
========================
ACMG/AMP Variant Pathogenicity Classifier using Tavtigian Bayesian framework.
Gene-specific calibration for TPP1.

Authors: Lena Traczuk, Dawid Fleischer
"""

import sys
import pandas as pd
import numpy as np

# Config and constants
GENE_PARAMS = {
    "inheritance": "autosomal_recessive",
    "pm2_threshold": 0.005,
    "bs1_threshold": 0.01,
    "ba1_threshold": 0.05,
    "critical_domains": [
        {"name": "propeptide_region", "aa_start": 20, "aa_end": 195},
        {"name": "catalytic_domain", "aa_start": 368, "aa_end": 563},
    ],
}

ACMG_POINTS = {
    "PVS1": 8, "PVS1_Strong": 4, "PVS1_Moderate": 2, "PVS1_Supporting": 1,
    "PS": 4, "PS3": 4, "PS1": 4,
    "PM": 2, "PM1": 2, "PM2": 1, "PM2_Supporting": 1, "PM4": 2, "PM5": 2,
    "PP": 1, "PP3": 1, "PP3_Strong": 2, "PP3_Moderate": 2, "PP3_Supporting": 1,
    "BP4_Supporting": -1, "BP4_Moderate": -2,
    "BP": -1, "BM": -2, "BS": -4, "BS1": -4, "BA1": -8, "BA": -8
}

ACMG_THRESHOLDS = {
    "pathogenic": 10,
    "likely_pathogenic": 6,
    "vus_upper": 5,
    "likely_benign_lower": -1,
    "benign": -6
}

def get_af(row):
    val = row.get("GNOMAD_Allele_Frequency", 0.0)
    if pd.isna(val) or val == "" or val == ".":
        return 0.0
    try:
        return float(val)
    except ValueError:
        return 0.0

def get_protein_position(row):
    val = row.get("Protein_position", "")
    if pd.isna(val) or val == "" or val == ".":
        return None
    try:
        if "-" in str(val):
            return int(str(val).split("-")[0])
        return int(val)
    except ValueError:
        return None

def is_in_critical_domain(pos, domains):
    if pos is None:
        return False
    for d in domains:
        if d["aa_start"] <= pos <= d["aa_end"]:
            return True
    return False

def classify_variant(row, params):
    criteria = []
    points = 0
    
    # 1. BA1 & BS1 (Population Frequencies)
    af = get_af(row)
    if af > params["ba1_threshold"]:
        criteria.append("BA1")
        points += ACMG_POINTS["BA1"]
    elif af > params["bs1_threshold"]:
        criteria.append("BS1")
        points += ACMG_POINTS["BS1"]
    
    # Skip pathogenic rules if BA1 is active (Stand-alone benign)
    if "BA1" in criteria:
        return {
            "acmg_criteria": ";".join(criteria),
            "acmg_points": points,
            "acmg_classification": "Benign"
        }

    # 2. PVS1 (Loss of Function)
    consequence = str(row.get("Consequence", ""))
    impact = str(row.get("IMPACT", ""))
    exon = str(row.get("EXON", ""))
    
    is_lof = any(c in consequence for d in ["nonsense", "frameshift", "splice_donor", "splice_acceptor"] for c in [d]) or impact == "HIGH"
    
    if is_lof:
        # Check if in last exon (NMD escape)
        try:
            if "/" in exon:
                current, total = map(int, exon.split("/"))
                if current == total:
                    criteria.append("PVS1_Strong")
                    points += ACMG_POINTS["PVS1_Strong"]
                else:
                    criteria.append("PVS1")
                    points += ACMG_POINTS["PVS1"]
            else:
                criteria.append("PVS1")
                points += ACMG_POINTS["PVS1"]
        except Exception:
            criteria.append("PVS1")
            points += ACMG_POINTS["PVS1"]

    # 3. PM2 (Rarity)
    if af == 0.0 or af < params["pm2_threshold"]:
        criteria.append("PM2_Supporting")
        points += ACMG_POINTS["PM2_Supporting"]

    # 4. PM1 (Critical Domain)
    prot_pos = get_protein_position(row)
    if is_in_critical_domain(prot_pos, params["critical_domains"]) and "missense" in consequence:
        criteria.append("PM1")
        points += ACMG_POINTS["PM1"]

    # 5. PP3 / BP4 (In silico predictors - REVEL & SpliceAI)
    revel = row.get("REVEL_score", ".")
    spliceai = row.get("SPLICEAI_DS_MAX", ".")
    
    revel_val = None
    if revel != "." and not pd.isna(revel):
        try:
            revel_val = float(revel)
        except ValueError:
            pass
            
    spliceai_val = None
    if spliceai != "." and not pd.isna(spliceai):
        try:
            spliceai_val = float(spliceai)
        except ValueError:
            pass

    # REVEL mapping
    if revel_val is not None:
        if revel_val >= 0.932:
            criteria.append("PP3_Strong")
            points += ACMG_POINTS["PP3_Strong"]
        elif revel_val >= 0.773:
            criteria.append("PP3_Moderate")
            points += ACMG_POINTS["PP3_Moderate"]
        elif revel_val >= 0.644:
            criteria.append("PP3_Supporting")
            points += ACMG_POINTS["PP3_Supporting"]
        elif revel_val <= 0.183:
            criteria.append("BP4_Moderate")
            points += ACMG_POINTS["BP4_Moderate"]
        elif revel_val <= 0.290:
            criteria.append("BP4_Supporting")
            points += ACMG_POINTS["BP4_Supporting"]

    # SpliceAI mapping (only if PP3 not already Strong/Moderate to avoid double-counting)
    if spliceai_val is not None and not any(c in ["PP3_Strong", "PP3_Moderate"] for c in criteria):
        if spliceai_val >= 0.8:
            criteria.append("PP3_Strong")
            points += ACMG_POINTS["PP3_Strong"]
        elif spliceai_val >= 0.5:
            criteria.append("PP3_Moderate")
            points += ACMG_POINTS["PP3_Moderate"]
        elif spliceai_val >= 0.2:
            criteria.append("PP3_Supporting")
            points += ACMG_POINTS["PP3_Supporting"]

    # 6. ClinVar proxies (PS3/PS1)
    clinvar_sig = str(row.get("CLINVAR_CLNSIG", "")).lower()
    if "pathogenic" in clinvar_sig and not "conflicting" in clinvar_sig:
        criteria.append("PS3")
        points += ACMG_POINTS["PS3"]
    elif "benign" in clinvar_sig and not "conflicting" in clinvar_sig:
        criteria.append("BS")
        points += ACMG_POINTS["BS"]

    # Final Classification
    if points >= ACMG_THRESHOLDS["pathogenic"]:
        classification = "Pathogenic"
    elif points >= ACMG_THRESHOLDS["likely_pathogenic"]:
        classification = "Likely_Pathogenic"
    elif points >= 0:
        classification = "VUS"
    elif points >= ACMG_THRESHOLDS["likely_benign_lower"]:
        classification = "Likely_Benign"
    else:
        classification = "Benign"

    return {
        "acmg_criteria": ";".join(criteria),
        "acmg_points": points,
        "acmg_classification": classification
    }

def main(in_path, out_path):
    df = pd.read_csv(in_path, sep="\t")
    results = []
    for _, row in df.iterrows():
        res = classify_variant(row, GENE_PARAMS)
        results.append(res)
    
    res_df = pd.DataFrame(results)
    final_df = pd.concat([df, res_df], axis=1)
    final_df.to_csv(out_path, sep="\t", index=False)
    print(f"Classified {len(final_df)} variants. Saved: {out_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: classify_acmg.py <input_tsv> <output_tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
