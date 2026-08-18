# omicsApp deployment image.
#
# Build from the repository root, not from this directory -- the COPY
# steps below need both packages in the build context:
#
#   docker build -t omicsapp:1.0 -f deploy/docker/Dockerfile .
#
# Layer order is chosen so that editing R code does not re-run the slow
# steps: system deps, then CRAN, then Bioconductor, then the gene-set
# prewarm, and only then the two packages, which change every commit.

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.source="https://github.com/puweilin/omicsApp"
LABEL org.opencontainers.image.description="Multi-omics analysis web app (omicsCore + omicsApp)"

# ---- 1. system dependencies -------------------------------------------
# pandoc renders the HTML report. The Bioconductor base image already
# carries the compilers and headers the Bioc stack needs.
RUN apt-get update && apt-get install -y --no-install-recommends \
        pandoc \
    && rm -rf /var/lib/apt/lists/*

# ---- 2. CRAN dependencies ---------------------------------------------
# Imports *and* Suggests. Suggests are not optional here: omicsCore
# gates features behind requireNamespace(), so a missing Suggest does
# not error, it silently removes the feature.
RUN R -q -e 'install.packages(c( \
      "shiny","bslib","bsicons","sass","htmltools","ggplot2","plotly","DT", \
      "shinyWidgets","shinyjs","promises","future", \
      "tibble","dplyr","stringr","readxl","openxlsx","jsonlite","scales", \
      "rlang","msigdbr","qs2", \
      "ggrepel","patchwork","ggpubr","circlize","here","knitr","rmarkdown", \
      "missForest","ActivePathways","pkgload"), \
      Ncpus = parallel::detectCores())'

# ---- 3. Bioconductor dependencies -------------------------------------
RUN R -q -e 'BiocManager::install(c( \
      "limma","clusterProfiler","DESeq2","edgeR","S4Vectors","enrichplot", \
      "fgsea","GSVA","ComplexHeatmap","impute","pcaMethods","tximport", \
      "GenomicFeatures"), ask = FALSE, update = FALSE, Ncpus = parallel::detectCores())'

# ---- 4. PDF reports (optional) ----------------------------------------
# Adds ~1.5 GB. Comment out if HTML reports are enough; the Report view
# gates the PDF button on rmarkdown being able to find a LaTeX engine.
RUN apt-get update && apt-get install -y --no-install-recommends \
        texlive-latex-base texlive-fonts-recommended texlive-latex-extra \
    && rm -rf /var/lib/apt/lists/*

# ---- 5. gene-set cache -------------------------------------------------
# The single most expensive runtime operation in the app is the first
# msigdbr() call (~10s). Baking the tables in makes it ~0.3s, for every
# user and every container.
ENV OMICSCORE_GENESET_CACHE=/opt/genesets
COPY deploy/docker/prewarm_genesets.R /tmp/prewarm_genesets.R
RUN Rscript /tmp/prewarm_genesets.R && rm /tmp/prewarm_genesets.R

# ---- 6. the packages ---------------------------------------------------
# omicsCore first: omicsApp Imports it. The `Remotes:` line in
# omicsApp/DESCRIPTION is only read by pak/devtools, so installing in
# order is all that is required here.
COPY packages/omicsCore /tmp/omicsCore
RUN R CMD INSTALL --no-docs /tmp/omicsCore && rm -rf /tmp/omicsCore
COPY packages/omicsApp /tmp/omicsApp
RUN R CMD INSTALL --no-docs /tmp/omicsApp && rm -rf /tmp/omicsApp

# ---- 7. runtime --------------------------------------------------------
# uid 1001 must match the owner of the per-user host directories that
# ShinyProxy bind-mounts at /data (see deploy/scripts/add_user.sh).
RUN useradd --create-home --uid 1001 omics \
    && mkdir -p /data && chown omics:omics /data

# One BLAS thread per worker. Without this each of the `workers`
# background processes opens as many threads as there are cores, and a
# handful of containers oversubscribe the host.
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OMICSAPP_DATA_DIR=/data \
    OMICSAPP_QUOTA_GB=50

USER omics
WORKDIR /home/omics
EXPOSE 3838

# host must be 0.0.0.0 for ShinyProxy to reach the container.
CMD ["R", "-q", "-e", "omicsApp::launch(host = '0.0.0.0', port = 3838, launch.browser = FALSE, workers = 2)"]
