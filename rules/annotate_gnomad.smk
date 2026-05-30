# rules/annotate_gnomad.smk

rule annotate_gnomad:
    input:
        vcf = "results/TPP1/TPP1_spliceai.vcf"
    output:
        tsv = "results/TPP1/TPP1_gnomad.tsv"
    conda:
        "../envs/python.yaml"
    shell:
        "python3 scripts/annotate_gnomad.py {input.vcf} {output.tsv}"
