"""
tests/test_pipeline.py
======================
Unit and integration tests for the TPP1 variant pathogenicity pipeline.

Authors: Lena Traczuk, Dawid Fleischer
"""

import pytest
import pandas as pd
import sys
from pathlib import Path

# Add scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from classify_acmg import classify_variant, get_af, get_protein_position, is_in_critical_domain

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

def test_get_af():
    row = pd.Series({"GNOMAD_Allele_Frequency": "0.001"})
    assert get_af(row) == 0.001
    
    row_empty = pd.Series({"GNOMAD_Allele_Frequency": "."})
    assert get_af(row_empty) == 0.0

def test_get_protein_position():
    row = pd.Series({"Protein_position": "208"})
    assert get_protein_position(row) == 208

def test_is_in_critical_domain():
    assert is_in_critical_domain(420, GENE_PARAMS["critical_domains"]) is True
    assert is_in_critical_domain(50, GENE_PARAMS["critical_domains"]) is True
    assert is_in_critical_domain(300, GENE_PARAMS["critical_domains"]) is False

def test_classify_variant():
    # Pathogenic variant
    row = pd.Series({
        "Consequence": "stop_gained",
        "IMPACT": "HIGH",
        "EXON": "6/13",
        "GNOMAD_Allele_Frequency": "0.0",
        "Protein_position": "208",
        "REVEL_score": ".",
        "SPLICEAI_DS_MAX": ".",
        "CLINVAR_CLNSIG": "Pathogenic",
    })
    res = classify_variant(row, GENE_PARAMS)
    assert res["acmg_classification"] in ["Pathogenic", "Likely_Pathogenic"]
