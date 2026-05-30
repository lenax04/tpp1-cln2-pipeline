# Snakefile
# =========
# Snakemake workflow for TPP1 variant classification.
# Authors: Lena Traczuk, Dawid Fleischer

configfile: "config/config.yaml"

rule all:
    input:
        config["classified_tsv"],
        "results/TPP1/figures/TPP1_classification_summary.png",
        "results/comparison/platform_comparison.png",
        "results/report/pipeline_report.html"

include: "rules/normalize.smk"
include: "rules/annotate_vep.smk"
include: "rules/annotate_spliceai.smk"
include: "rules/annotate_gnomad.smk"
include: "rules/annotate_clinvar.smk"
include: "rules/classify_acmg.smk"
include: "rules/compare_platforms.smk"
include: "rules/report.smk"
