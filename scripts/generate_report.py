#!/usr/bin/env python3
"""
scripts/generate_report.py
==========================
Generate HTML summary report for the TPP1 variant classification pipeline.

Authors: Lena Traczuk, Dawid Fleischer
"""

import sys
import pandas as pd
from pathlib import Path

def main(in_path, comp_path, out_html):
    df = pd.read_csv(in_path, sep="\t")
    comp_df = pd.read_csv(comp_path, sep="\t")
    
    # Calculate statistics
    total = len(df)
    counts = df["acmg_classification"].value_counts()
    p_lp = counts.get("Pathogenic", 0) + counts.get("Likely_Pathogenic", 0)
    vus = counts.get("VUS", 0)
    b_lb = counts.get("Benign", 0) + counts.get("Likely_Benign", 0)
    
    # HTML content
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>MicroSnake Variant Classification Report — TPP1</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 30px; background-color: #f8f9fa; color: #333; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }}
        h2 {{ color: #16a085; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 15px; background: white; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #2c3e50; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .card {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 20px; }}
        .flex {{ display: flex; gap: 20px; }}
        .stat {{ flex: 1; text-align: center; padding: 15px; border-radius: 6px; color: white; font-weight: bold; }}
        .stat-total {{ background-color: #34495e; }}
        .stat-plp {{ background-color: #e74c3c; }}
        .stat-vus {{ background-color: #f1c40f; }}
        .stat-blb {{ background-color: #2ecc71; }}
    </style>
</head>
<body>
    <h1>🧬 MicroSnake Variant Classification Report — TPP1 (CLN2)</h1>
    <p><strong>Authors:</strong> Lena Traczuk & Dawid Fleischer</p>
    <p><strong>Date:</strong> 2026-05-30</p>
    
    <div class="card flex">
        <div class="stat stat-total">
            <h3>Total Variants</h3>
            <p style="font-size: 24px; margin: 5px 0;">{total}</p>
        </div>
        <div class="stat stat-plp">
            <h3>Pathogenic / Likely Pathogenic</h3>
            <p style="font-size: 24px; margin: 5px 0;">{p_lp} ({p_lp/total*100:.1f}%)</p>
        </div>
        <div class="stat stat-vus">
            <h3>VUS</h3>
            <p style="font-size: 24px; margin: 5px 0;">{vus} ({vus/total*100:.1f}%)</p>
        </div>
        <div class="stat stat-blb">
            <h3>Benign / Likely Benign</h3>
            <p style="font-size: 24px; margin: 5px 0;">{b_lb} ({b_lb/total*100:.1f}%)</p>
        </div>
    </div>

    <div class="card">
        <h2>📊 Classification Distribution</h2>
        <table>
            <tr>
                <th>ACMG Classification</th>
                <th>Count</th>
                <th>Percentage</th>
            </tr>
    """
    
    for cls in ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]:
        cnt = counts.get(cls, 0)
        html += f"""
            <tr>
                <td><strong>{cls}</strong></td>
                <td>{cnt}</td>
                <td>{cnt/total*100:.1f}%</td>
            </tr>
        """
        
    html += f"""
        </table>
    </div>

    <div class="card">
        <h2>🔍 Platform Concordance (vs. Franklin by Genoox)</h2>
        <table>
            <tr>
                <th>Gene</th>
                <th>Disease</th>
                <th>Total Variants</th>
                <th>Pipeline P/LP</th>
                <th>Concordant P/LP</th>
                <th>Pipeline Unique</th>
                <th>Franklin Unique</th>
                <th>VUS Elevated by Franklin</th>
                <th>Concordance (%)</th>
            </tr>
            <tr>
                <td>{comp_df.iloc[0]['Gene']}</td>
                <td>{comp_df.iloc[0]['Disease']}</td>
                <td>{comp_df.iloc[0]['Total_Variants']}</td>
                <td>{comp_df.iloc[0]['Pipeline_P_LP_Total']}</td>
                <td>{comp_df.iloc[0]['Concordant_P_LP']}</td>
                <td>{comp_df.iloc[0]['Pipeline_Unique_P_LP']}</td>
                <td>{comp_df.iloc[0]['Franklin_Unique_P_LP']}</td>
                <td>{comp_df.iloc[0]['VUS_Elevated_by_Franklin']}</td>
                <td><strong>{comp_df.iloc[0]['Concordance_Pct']}%</strong></td>
            </tr>
        </table>
    </div>
</body>
</html>
    """
    
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    with open(out_html, "w") as f:
        f.write(html)
    print(f"Generated HTML report: {out_html}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: generate_report.py <input_tsv> <comp_tsv> <output_html>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
