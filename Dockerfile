FROM condaforge/mambaforge:23.3.1-1

LABEL maintainer="Lena Traczuk <lena.traczuk@example.com>, Dawid Fleischer <dawid.fleischer@example.com>"
LABEL description="TPP1/CLN2 Variant Pathogenicity Classification Pipeline"
LABEL version="1.0.0"

RUN apt-get update && apt-get install -y wget curl git tabix bcftools && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN mamba install -y -c conda-forge -c bioconda snakemake=7.32 python=3.11 pandas>=2.0 numpy>=1.24 matplotlib>=3.7 seaborn>=0.12 scipy>=1.10 pytest && mamba clean --all -y

WORKDIR /pipeline
COPY . .
CMD ["snakemake", "--dry-run", "--cores", "1"]
