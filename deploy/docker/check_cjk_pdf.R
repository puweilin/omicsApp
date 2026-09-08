# Does this image actually render a Chinese project name into a PDF?
#
# Asked at build time because the failure mode is quiet. Without the font
# packages XeLaTeX does not stop -- it typesets empty boxes -- so the first
# person to find out would be a user downloading a report of blank
# rectangles. And without xeCJK the engine falls back to a font that has no
# CJK coverage, with the same result.
#
# The old pdfLaTeX default failed loudly instead, with "Unicode character
# 蛋 (U+86CB) not set up for use with LaTeX", which is at least honest. The
# point of this check is to keep that honesty after switching engines.
#
# ---- how the embedding is verified -------------------------------------
#
# With pdffonts, not by searching the file for the font's name. XeLaTeX
# writes its cross-reference table and font dictionaries into compressed
# object streams, so a PDF it produced contains neither "Noto" nor even
# the literal "BaseFont" as plain bytes -- a raw grep reports every such
# PDF as having no embedded font, whatever is really in it. That check
# could only ever fail, which is worse than no check: it fails the build
# on a correct image and teaches whoever sees it to stop reading.
#
# pdffonts reads the structure and prints one row per font with an `emb`
# column. That is the question being asked, answered by the tool that
# exists for it.

suppressPackageStartupMessages(library(rmarkdown))

if (!rmarkdown::pandoc_available()) {
  stop("pandoc is not available; the report view cannot render anything.")
}
if (!nzchar(Sys.which("pdffonts"))) {
  stop("pdffonts is not available; install poppler-utils. ",
       "The build cannot verify that the CJK font was embedded without it.")
}

dir <- file.path(tempdir(), "cjk-check")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
rmd <- file.path(dir, "cjk.Rmd")

writeLines(c(
  "---",
  'title: "蛋白组项目 · Proteomics report"',
  "output:",
  "  pdf_document:",
  "    latex_engine: xelatex",
  "header-includes:",
  "  - \\usepackage{xeCJK}",
  "  - \\setCJKmainfont{Noto Sans CJK SC}",
  # Named so xeCJK does not warn about an undefined mono family. The
  # warning is harmless and the noise is not: a build log that always
  # carries warnings is one nobody reads when a real one appears.
  "  - \\setCJKmonofont{Noto Sans Mono CJK SC}",
  "---",
  "",
  "样本数量与差异蛋白：这一行必须能排版出来。",
  "",
  "```{r}",
  "knitr::kable(data.frame(样本 = c('对照组', '实验组'), n = c(12L, 12L)))",
  "```"
), rmd)

out <- rmarkdown::render(rmd, output_file = "cjk.pdf", output_dir = dir,
                         quiet = TRUE)

size <- file.info(out)$size
if (!file.exists(out) || is.na(size) || size < 1000) {
  stop("CJK PDF check produced no usable file.")
}

# A PDF that rendered blank boxes is the same size as one that rendered
# glyphs, so size is not the test.
fonts <- system2("pdffonts", shQuote(out), stdout = TRUE, stderr = TRUE)

# Columns are name, type, encoding, emb, sub, uni, object ID. Subset
# fonts carry a six-letter tag -- `ABCDEF+NotoSansCJKsc-Regular` -- so the
# name is matched rather than compared.
rows <- grep("Noto|CJK", fonts, value = TRUE, ignore.case = TRUE)
embedded <- grep("\\byes\\b", rows, value = TRUE)

if (length(embedded) == 0L) {
  stop("No CJK font is embedded; the report would print empty boxes. ",
       "Check that fonts-noto-cjk and texlive-lang-chinese are installed.\n",
       "pdffonts said:\n", paste(fonts, collapse = "\n"))
}

cat("CJK PDF check passed:", size, "bytes\n")
cat("  embedded:", sub("\\s.*$", "", trimws(embedded[[1L]])), "\n")
