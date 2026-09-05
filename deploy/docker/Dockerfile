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
#
# `install.packages()` warns on failure rather than erroring, so without
# the check below a package that will not build produces a *successful*
# layer and an image quietly missing a feature. Fail here instead, where
# the log still says which package and why.
#
# ggpubr is deliberately absent. It pulls Deriv -> doBy -> pbkrtest ->
# car -> rstatix, and Deriv's current source needs R_ClosureFormals(),
# an R 4.5 C API this image's R 4.4.2 does not have. Nothing calls
# ggpubr: plots use theme_analysis()/theme_volcano() from the framework.
#
# CRAN is pinned, because Bioconductor 3.20 is archived and today's CRAN
# has moved nearly two years past it. Left on `latest` the mismatch
# surfaces one package at a time: Deriv wants R_ClosureFormals(), an R
# 4.5 C API this R 4.4.2 lacks, and ggtree 3.14.0 calls
# getFromNamespace("check_linewidth", "ggplot2"), which ggplot2 4.0
# removed. Pinning fixes the class rather than the instances, and makes
# the build reproducible as a side effect.
#
# The date is from the *end* of 3.20's life, not its release. The mirror
# serves the final state of the 3.20 branch, which kept moving all cycle:
# a snapshot from two days after 3.20 shipped had all 228 packages
# present and still failed, because BiocParallel by then wanted
# BH >= 1.87.0, released later. 2025-03-01 is the latest date that
# satisfies every requirement while msigdbr is still < 10 -- see WATCH in
# deploy/scripts/check_pins.R, which verifies both and should be run
# before this line is changed.
ARG CRAN_SNAPSHOT=2025-03-01
ARG BIOC_MIRROR=https://mirrors.westlake.edu.cn/bioconductor
ENV CRAN=https://p3m.dev/cran/__linux__/noble/${CRAN_SNAPSHOT}
COPY deploy/docker/setup_repos.R /tmp/setup_repos.R
RUN Rscript /tmp/setup_repos.R && rm /tmp/setup_repos.R

RUN R -q -e 'pkgs <- c( \
      "shiny","bslib","bsicons","sass","htmltools","ggplot2","plotly","DT", \
      "shinyWidgets","shinyjs","promises","future", \
      "tibble","dplyr","stringr","readxl","openxlsx","jsonlite","scales", \
      "rlang","msigdbr","qs2", \
      "ggrepel","patchwork","circlize","here","knitr","rmarkdown", \
      "imputeLCMD","ActivePathways","pkgload"); \
      install.packages(pkgs, Ncpus = parallel::detectCores()); \
      missing <- setdiff(pkgs, rownames(installed.packages())); \
      if (length(missing)) stop("CRAN packages failed to install: ", \
                                paste(missing, collapse = ", "))'

# ---- 3. Bioconductor dependencies -------------------------------------
# `vsn` is not optional in practice: it is the default normalisation in
# the import view and the method the legacy proteomics pipeline has
# always used. Without it every proteomics import fails on the default
# setting with "Package 'vsn' is required".
#
# Checked for the same reason as the CRAN layer above, and it matters
# more here: none of these have binaries on p3m, so every one is a
# source build with its own way to fail.
#
# The repositories are named explicitly rather than left to
# BiocManager::install(), because 3.20 is now an archived release and
# bioconductor.org redirects it to an OSN bucket at XSEDE that is slow
# to unreachable. The failure mode is the bad one: no error, just a
# download waiting on a host that never answers, so the build appears to
# hang after printing "Replacement repositories".
#
# data/annotation is in the list because clusterProfiler needs GO.db.
# Override BIOC_MIRROR to build somewhere the upstream is fast:
#   docker build --build-arg BIOC_MIRROR=https://bioconductor.org ...
#
# 12 packages pull in 168: clusterProfiler alone accounts for 120,
# enrichplot 118 and GSVA 101, overlapping heavily. That closure, all
# built from source because p3m carries no Bioconductor binaries, is
# what the 30-45 minutes here is.
#
# tximport and GenomicFeatures are deliberately absent. Neither is ever
# called: diff-deseq2.R and diff-edger.R read a length matrix the caller
# may have left at input$misc$tximport$length, which is metadata, not a
# package call, and there is no Salmon import path yet. Together they
# cost 11 further packages of which rtracklayer, Rsamtools, Rhtslib,
# GenomicAlignments and RCurl are all C/C++ builds. They stay in
# omicsCore's Suggests and install-optional.R's rnaseq bundle, for
# whoever installs the package outside this image.
# Split into three, ordered by how deep their dependency chains run, so
# that a failure in one keeps the earlier ones cached. That is worth the
# extra layers: this closure is 109 source builds, both failures so far
# landed in the last two groups -- BiocParallel under DESeq2, ggtree
# under clusterProfiler -- and each time the whole layer was discarded,
# because Docker commits only successful ones.
#
# Repositories come from Rprofile.site, written above.

# 3a. Foundations: 24 packages, none deep.
RUN R -q -e 'pkgs <- c("S4Vectors","limma","impute","pcaMethods","edgeR", \
                       "vsn","ComplexHeatmap"); \
      install.packages(pkgs, Ncpus = parallel::detectCores()); \
      missing <- setdiff(pkgs, rownames(installed.packages())); \
      if (length(missing)) stop("Bioconductor 3a failed: ", \
                                paste(missing, collapse = ", "))'

# 3b. 23 more, and the first heavy C++ builds (BiocParallel, Rhdf5lib).
RUN R -q -e 'pkgs <- c("fgsea","DESeq2"); \
      install.packages(pkgs, Ncpus = parallel::detectCores()); \
      missing <- setdiff(pkgs, rownames(installed.packages())); \
      if (length(missing)) stop("Bioconductor 3b failed: ", \
                                paste(missing, collapse = ", "))'

# 3c. The long tail: 62 packages, most of the time, and the whole
# clusterProfiler/enrichplot/ggtree chain.
RUN R -q -e 'pkgs <- c("GSVA","enrichplot","clusterProfiler"); \
      install.packages(pkgs, Ncpus = parallel::detectCores()); \
      missing <- setdiff(pkgs, rownames(installed.packages())); \
      if (length(missing)) stop("Bioconductor 3c failed: ", \
                                paste(missing, collapse = ", "))'

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
# The container writes into a host directory through a bind mount, and
# the kernel decides whether it may by comparing *numbers*: the uid of
# the process against the uid owning the directory. Names do not enter
# into it -- "omics" exists only inside this image.
#
# With root on the host you leave this alone and chown the directories
# to 1001 (deploy/scripts/add_user.sh). Without root you cannot chown to
# a uid you do not own, so build the image as your own uid instead and
# the numbers match without any chown:
#
#   docker build --build-arg APP_UID=$(id -u) -t omicsapp:1.0 -f deploy/docker/Dockerfile .
#
# Getting this wrong is quiet: the app reports "could not save" and
# carries on, so the user's work simply never persists.
ARG APP_UID=1001
RUN useradd --create-home --uid ${APP_UID} omics \
    && mkdir -p /data /data/raw \
    && chown -R omics:omics /data

# One BLAS thread per worker, and no nested parallelism.
#
# ENV alone does not achieve this, which is worth stating plainly
# because the previous version of this file claimed it did: R reads
# etc/Renviron.site during startup and its assignments *override* the
# inherited environment. The base image sets OMP_NUM_THREADS=2,
# OMP_THREAD_LIMIT=2 and BIOCPARALLEL_WORKER_NUMBER=4 there, so every R
# process got those regardless of what ENV said. Appending wins, because
# R applies the file in order.
#
# The base image's own comments say those values work around package
# check timeouts -- "a BayesSpace TIMEOUT on Linux", "a snifter TIMEOUT
# on Mac". None of that applies to an application container that builds
# nothing and runs no checks.
#
# What does apply is that launch() uses future::multicore, which forks.
# BiocParallel's default back-end forks too, so DESeq2, GSVA and fgsea
# would fork again inside an already-forked worker: 2 workers x 4 = 8
# processes on a 2-CPU container. Worse, OpenMP and fork() are a known
# hazard -- a process that has started an OpenMP thread pool and then
# forks leaves the child with a pool in an undefined state, and the
# symptom is a session that hangs with nothing on screen. One thread in
# the parent is the standard way to avoid it.
#
# The cost is that a lone analysis is single-threaded. That is the right
# side of the trade here: the container has 2 CPUs, run_diff() takes
# seconds, and a hang has no error message to debug.
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OMICSAPP_DATA_DIR=/data \
    OMICSAPP_QUOTA_GB=50

# A UTF-8 locale, stated rather than inherited. The rocker base sets
# one, but nothing here depended on that being true until it was tested:
# under LC_ALL=C, R cannot hand a path with a non-ASCII character to the
# operating system, so a project named in Chinese could be neither saved
# nor opened, and a sample name from a CSV reached the browser as bytes.
# The code now survives a C locale (see the c-locale tests in both
# packages); the image still should not run in one.
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8
RUN printf '%s\n' \
      '' \
      '# Appended by deploy/docker/Dockerfile; overrides the base image.' \
      'OMP_NUM_THREADS=1' \
      'OMP_THREAD_LIMIT=1' \
      'KMP_DEVICE_THREAD_LIMIT=1' \
      'KMP_TEAMS_THREAD_LIMIT=1' \
      'BIOCPARALLEL_WORKER_NUMBER=1' \
    >> /usr/local/lib/R/etc/Renviron.site \
 && R -q --no-echo -e 'v <- c("OMP_NUM_THREADS","OMP_THREAD_LIMIT","BIOCPARALLEL_WORKER_NUMBER"); \
      got <- Sys.getenv(v); \
      cat(paste0("  ", v, " = ", got), sep = "\n"); \
      if (!all(got == "1")) stop("Renviron.site override did not take effect")'

USER omics
WORKDIR /home/omics
EXPOSE 3838

# host must be 0.0.0.0 for ShinyProxy to reach the container.
CMD ["R", "-q", "-e", "omicsApp::launch(host = '0.0.0.0', port = 3838, launch.browser = FALSE, workers = 2)"]
