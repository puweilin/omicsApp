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

suppressPackageStartupMessages(library(rmarkdown))

if (!rmarkdown::pandoc_available()) {
  stop("pandoc is not available; the report view cannot render anything.")
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
# glyphs, so size is not the test. The font has to be embedded, and its
# name appears in the file when it is.
raw <- readBin(out, "raw", n = size)
txt <- rawToChar(raw[raw != as.raw(0)])
if (!grepl("Noto|CJK", txt, useBytes = TRUE)) {
  stop("The CJK font was not embedded; the report would print empty boxes. ",
       "Check that fonts-noto-cjk and texlive-lang-chinese are installed.")
}

cat("CJK PDF check passed:", size, "bytes\n")
