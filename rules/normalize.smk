# rules/normalize.smk

rule normalize_variants:
    input:
        vcf = config["raw_vcf"]
    output:
        vcf = config["normalized_vcf"]
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools norm -m -any {input.vcf} | "
        "bcftools norm -f GRCh38.fa -o {output.vcf}"
