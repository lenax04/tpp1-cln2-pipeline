# rules/report.smk

rule generate_report:
    input:
        tsv = config["classified_tsv"],
        comp = "results/comparison/concordance_summary.tsv"
    output:
        html = config["report_html"]
    conda:
        "../envs/python.yaml"
    shell:
        "python3 scripts/generate_report.py {input.tsv} {input.comp} {output.html}"
