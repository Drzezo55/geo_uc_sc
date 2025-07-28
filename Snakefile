configfile: "config/config.yaml"

rule all:
    input:
        expand("output/{file}", file=config["output_files"])

rule run_seurat_pipeline:
    input:
        script="scripts/script.R"
    output:
        expand("output/{file}", file=config["output_files"])
    conda:
        "envs/snakemake_env.yaml" # <--- Add this line, pointing to your environment file
    shell:
        "Rscript {input.script}"
