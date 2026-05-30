# rules/annotate_spliceai.smk

rule annotate_spliceai:
    input:
        vcf = config["annotated_vcf"]
    output:
        vcf = "results/TPP1/TPP1_spliceai.vcf"
    conda:
        "../envs/spliceai.yaml"
    shell:
        "spliceai -I {input.vcf} -O {output.vcf} -R GRCh38.fa -D 500"
