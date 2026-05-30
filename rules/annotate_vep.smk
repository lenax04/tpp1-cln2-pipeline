# rules/annotate_vep.smk

rule annotate_vep:
    input:
        vcf = config["normalized_vcf"]
    output:
        vcf = config["annotated_vcf"]
    conda:
        "../envs/vep.yaml"
    shell:
        "vep --input_file {input.vcf} --output_file {output.vcf} --vcf "
        "--cache --dir_cache local_cache --assembly GRCh38 --pick "
        "--plugin dbNSFP,REVEL,BayesDel_noAF --check_existing"
