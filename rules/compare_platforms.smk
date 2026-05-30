# rules/compare_platforms.smk

rule compare_platforms:
    input:
        tsv = config["classified_tsv"]
    output:
        fig = "results/comparison/platform_comparison.png",
        tsv = "results/comparison/concordance_summary.tsv"
    conda:
        "../envs/python.yaml"
    shell:
        "python3 scripts/compare_platforms.py {input.tsv} {output.fig} {output.tsv}"
