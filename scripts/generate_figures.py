#!/usr/bin/env python3
"""
scripts/generate_figures.py
===========================
Generate publication-quality figures for TPP1 variant classification.

Authors: Lena Traczuk, Dawid Fleischer
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def main(in_path, out_dir):
    df = pd.read_csv(in_path, sep="\t")
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    sns.set_theme(style="whitegrid")
    
    # Figure 1: Classification Summary
    plt.figure(figsize=(8, 6))
    order = ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]
    colors = ["#d62728", "#ff7f0e", "#bcbd22", "#2ca02c", "#1f77b4"]
    
    ax = sns.countplot(
        data=df,
        x="acmg_classification",
        order=order,
        palette=colors
    )
    plt.title("Distribution of ACMG/AMP Variant Classifications for TPP1", fontsize=14, fontweight="bold")
    plt.xlabel("Classification", fontsize=12)
    plt.ylabel("Variant Count", fontsize=12)
    
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(
                f'{int(height)}',
                (p.get_x() + p.get_width() / 2., height),
                ha='center', va='bottom',
                fontsize=11, color='black',
                xytext=(0, 3),
                textcoords='offset points'
            )
    plt.tight_layout()
    plt.savefig(out_path / "TPP1_classification_summary.png", dpi=300)
    plt.close()
    
    # Figure 2: Variant Spectrum (Consequence vs Classification)
    plt.figure(figsize=(10, 6))
    sns.countplot(
        data=df,
        y="Consequence",
        hue="acmg_classification",
        hue_order=order,
        palette=colors
    )
    plt.title("Variant Spectrum: Functional Consequences vs. ACMG Classifications", fontsize=14, fontweight="bold")
    plt.xlabel("Variant Count", fontsize=12)
    plt.ylabel("Molecular Consequence", fontsize=12)
    plt.legend(title="ACMG Class", loc="lower right")
    plt.tight_layout()
    plt.savefig(out_path / "TPP1_variant_spectrum.png", dpi=300)
    plt.close()
    
    # Figure 3: gnomAD AF distribution
    plt.figure(figsize=(8, 6))
    # Filter out 0 AF for log scale
    af_df = df[df["GNOMAD_Allele_Frequency"] > 0].copy()
    
    sns.stripplot(
        data=af_df,
        x="acmg_classification",
        y="GNOMAD_Allele_Frequency",
        order=order,
        hue="acmg_classification",
        palette=colors,
        jitter=0.25,
        size=6,
        alpha=0.7,
        legend=False
    )
    plt.yscale("log")
    plt.title("gnomAD Allele Frequency Distribution by ACMG Classification", fontsize=14, fontweight="bold")
    plt.xlabel("ACMG Classification", fontsize=12)
    plt.ylabel("gnomAD Allele Frequency (Log Scale)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path / "TPP1_gnomad_af_distribution.png", dpi=300)
    plt.close()
    
    print(f"Generated TPP1 figures in {out_dir}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: generate_figures.py <input_tsv> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
