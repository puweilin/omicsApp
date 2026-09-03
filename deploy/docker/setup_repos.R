#!/usr/bin/env Rscript
#
# Bake the pinned repositories into Rprofile.site, so every later layer
# and the running container agree on where packages come from.
#
# Reads CRAN and BIOC_MIRROR from the environment -- Docker exposes both
# ARG and ENV to a RUN -- and writes literal URLs, not another
# Sys.getenv(). ARG values do not survive into the running image, so a
# deferred lookup would leave the container with an empty repository.
#
# Doing this once rather than repeating options(repos=) in each install
# layer is what lets those layers be split: a failure in the last of
# them keeps the earlier ones cached, which matters because the
# Bioconductor closure is 109 source builds and a failure used to cost
# all of them.

bioc_mirror <- Sys.getenv("BIOC_MIRROR")
cran        <- Sys.getenv("CRAN")

if (!nzchar(bioc_mirror)) stop("BIOC_MIRROR is not set")
if (!nzchar(cran))        stop("CRAN is not set")

bioc <- paste0(bioc_mirror, "/packages/", Sys.getenv("BIOC_RELEASE", "3.20"))

lines <- c(
  "",
  "# Written by deploy/docker/setup_repos.R at image build time.",
  "options(",
  "  repos = c(",
  sprintf('    BioCsoft = "%s/bioc",', bioc),
  sprintf('    BioCann  = "%s/data/annotation",', bioc),
  sprintf('    BioCexp  = "%s/data/experiment",', bioc),
  sprintf('    CRAN     = "%s"', cran),
  "  ),",
  "  # Source tarballs over a mirror; R's 60s default is not generous.",
  "  timeout = 300",
  ")"
)

# R reads R_PROFILE in preference to etc/Rprofile.site, which is what
# makes this script testable outside a container: set it, run, and the
# fresh session below reads the same file the write went to.
path <- Sys.getenv("R_PROFILE", unset = file.path(R.home("etc"), "Rprofile.site"))
cat(lines, sep = "\n", file = path, append = TRUE)
cat("\n", file = path, append = TRUE)

# Prove it took effect in a *fresh* session rather than trusting the
# write: a typo here would otherwise surface as a confusing download
# failure three layers later.
out <- system2("R", c("-q", "--no-echo", "-e", shQuote(
  'r <- getOption("repos"); cat(paste(names(r), r, sep = " = "), sep = "\n")')),
  env = paste0("R_PROFILE=", path), stdout = TRUE)
cat("Repositories now in force:\n")
cat(paste0("  ", out), sep = "\n")

for (want in c("BioCsoft", "BioCann", "BioCexp", "CRAN")) {
  if (!any(grepl(paste0("^", want, " = http"), out))) {
    stop("Rprofile.site did not take effect: no ", want)
  }
}
