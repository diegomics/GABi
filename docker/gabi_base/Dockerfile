# minimal Debian base
FROM debian:trixie-slim

# build-time only (ARG does not persist into the running container)
ARG DEBIAN_FRONTEND=noninteractive
ARG MICROMAMBA_VERSION=2.9.0-0
# PLINK 2.0 Alpha 7.8 (AVX2 build, requires a CPU with AVX2)
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260919.zip
ARG IMAGE_VERSION=0.2.2

LABEL org.opencontainers.image.title="gabi_base" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.description="Base image for the GABi Snakemake pipeline"

# environment lives in /opt/conda; micromamba's package cache in /opt/micromamba
ENV CONDA_DIR=/opt/conda
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH="$CONDA_DIR/bin:$PATH"

# fail a RUN if any command in a pipe fails, not just the last one
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# 1. System dependencies and locales
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates curl unzip procps locales && \
    sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en

# 2. micromamba as a standalone static binary
RUN curl -fsSL -o /usr/local/bin/micromamba \
        "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" && \
    chmod +x /usr/local/bin/micromamba && \
    micromamba --version

# 3. Conda environment: single solve from the pinned spec
#    Flexible channel priority is required: the solution mixes bioconda and
#    conda-forge builds (e.g. bioconda perl modules on conda-forge deps).
COPY gabi_base_docker_env.yml /opt/gabi_base_docker_env.yml
RUN micromamba config set channel_priority flexible && \
    micromamba create -y -p "$CONDA_DIR" -f /opt/gabi_base_docker_env.yml && \
    micromamba env export -p "$CONDA_DIR" --explicit > /opt/gabi_base_docker_env.lock.txt && \
    micromamba clean -a -y

# 4. PLINK 2 standalone binary
RUN curl -fsSL -o /tmp/plink2.zip "$PLINK2_URL" && \
    unzip -o -q /tmp/plink2.zip -d /usr/local/bin/ && \
    rm /tmp/plink2.zip && \
    chmod +x /usr/local/bin/plink2

# 5. Verify: every executable exists, R packages load, Python modules import.
#    Reports all missing tools before failing.
RUN missing=0; \
    for t in micromamba python Rscript fastp fastqc NanoPlot multiqc seqtk seqkit \
             csvtk pigz bbduk.sh bbmap.sh gfastats bwa-mem2 minibwa minimap2 \
             samtools bamtools sambamba mosdepth bedops bedtools bcftools psmc \
             msmc2_Linux blastn makeblastdb esearch efetch enaDataGet.py enaGroupGet.py aria2c plink2; do \
        command -v "$t" >/dev/null || { echo "MISSING: $t"; missing=1; }; \
    done; \
    [ "$missing" -eq 0 ] && \
    Rscript -e 'for (p in c("RColorBrewer","ggplot2","dplyr","tidyr","scales","data.table","gridExtra","optparse")) suppressPackageStartupMessages(library(p, character.only = TRUE))' && \
    python -c "import numpy, scipy, tqdm, kaleido, plotly, pysam" && \
    echo "All tools present."

# Informational only: print versions to the build log (never fails the build)
RUN set +o pipefail; \
    for cmd in "fastp --version" "fastqc --version" "NanoPlot --version" \
               "multiqc --version" "samtools --version" "bcftools --version" \
               "sambamba --version" "bwa-mem2 version" "minimap2 --version" \
               "bedtools --version" "plink2 --version"; do \
        echo "--- $cmd: $($cmd 2>&1 | head -1)"; \
    done

CMD ["/bin/bash"]
