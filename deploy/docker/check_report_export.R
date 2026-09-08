# Does the Report view's download actually produce a file in this image?
#
# check_cjk_pdf.R already proved that XeLaTeX and the CJK fonts are here,
# and it passed while the application's own PDF button failed on the
# server. It rendered an Rmd it wrote into a temporary directory, as root,
# to a path ending in `.pdf`. Every one of those differs from what happens
# when a user clicks the button:
#
#   * the input is the template inside the installed package, under
#     /usr/local/lib/R/site-library -- owned by root, read-only to the
#     account this container runs as. LaTeX writes its .aux and .log
#     beside its input, so the render stopped with `I can't write on
#     file`. On a development checkout that directory is writable and
#     nothing was ever noticed.
#   * the process is `omics`, not root.
#   * Shiny's downloadHandler passes `content()` a bare tempfile(): no
#     extension, and inside tempdir() under the name `file<hex>`.
#
# So this check calls the real export_report(), as the real user, with the
# real path shape, after the packages are installed. It runs after
# `USER omics` for that reason -- a RUN before it would be root again and
# would prove nothing.

suppressPackageStartupMessages(library(omicsCore))

if (!nzchar(Sys.which("pdftotext"))) {
  stop("pdftotext is not available; install poppler-utils. The build ",
       "cannot read the report back to check it without it.")
}

set.seed(1)
n_feat <- 20L
n_samp <- 6L
expr <- matrix(
  rnorm(n_feat * n_samp, mean = 5, sd = 1),
  nrow = n_feat,
  dimnames = list(paste0("f", seq_len(n_feat)), paste0("s", seq_len(n_samp)))
)
meta <- data.frame(group = rep(c("ctrl", "case"), each = 3L),
                   row.names = paste0("s", seq_len(n_samp)))
feat <- data.frame(feature_id = paste0("f", seq_len(n_feat)),
                   feature_symbol = paste0("GENE", seq_len(n_feat)),
                   row.names = paste0("f", seq_len(n_feat)))

input <- omics_input(expr, meta, feat, omics_type = "proteomics",
                     assay_type = "normalized_intensity")
# Named as the projects here are named. A Latin-only name would not
# exercise xeCJK at all.
project <- omics_project("蛋白组项目",
                         experiments = list(proteo = input))

tpl_dir <- dirname(system.file("rmd", "default-report.Rmd",
                               package = "omicsCore"))
before <- list.files(tpl_dir, all.files = TRUE, no.. = TRUE)

for (fmt in c("html", "pdf")) {
  out <- tempfile()                     # no extension, exactly like Shiny
  export_report(project, out, format = fmt, overwrite = TRUE)

  if (!file.exists(out)) {
    stop("export_report(format = \"", fmt, "\") reported success but ",
         "wrote nothing to the path it was given. This is what the ",
         "download buttons serve, so they would serve an empty file.")
  }
  size <- file.info(out)$size
  if (is.na(size) || size < 1000) {
    stop("The ", fmt, " report is ", size, " bytes, which is not a report.")
  }
  cat(sprintf("  %-4s %d bytes\n", fmt, size))

  if (identical(fmt, "pdf")) {
    # Read the name back out of the finished PDF rather than looking for a
    # font by name in it. Which font is embedded is not the question --
    # `Noto Sans CJK SC` appears in the font table as `NotoSansCJKsc`, and
    # `Songti SC` as `STSongti-SC-Regular`, so a name match encodes a guess
    # about which font the image happens to ship and rejects images that
    # are perfectly fine.
    #
    # Whether the characters come back out is the question, and it catches
    # more than a missing font. Under a non-UTF-8 locale the name is
    # mangled long before LaTeX sees it -- correct fonts, embedded, and the
    # title still wrong -- which no font table would show.
    got <- system2("pdftotext", c(shQuote(out), "-"), stdout = TRUE,
                   stderr = FALSE)
    if (!any(grepl(project$name, got, fixed = TRUE))) {
      stop("The project name does not come back out of the exported PDF, ",
           "so the report would show empty boxes or mojibake where the ",
           "name should be.\nLooked for: ", project$name,
           "\nThe PDF's text begins:\n",
           paste(utils::head(got, 8L), collapse = "\n"),
           "\npdffonts:\n",
           paste(system2("pdffonts", shQuote(out), stdout = TRUE,
                         stderr = TRUE), collapse = "\n"))
    }
  }
  unlink(out)
}

after <- list.files(tpl_dir, all.files = TRUE, no.. = TRUE)
stray <- setdiff(after, before)
if (length(stray) > 0L) {
  stop("The render wrote into the installed package directory: ",
       paste(stray, collapse = ", "),
       ". That directory is read-only for some deployments, and this is ",
       "how the PDF export broke before.")
}

cat("Report export check passed (as ", Sys.info()[["user"]], ")\n", sep = "")
