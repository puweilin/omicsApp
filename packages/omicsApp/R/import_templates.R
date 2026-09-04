# Downloadable import templates.
#
# The import classifier is a heuristic over whatever a vendor happened to
# send. It is good at counts-and-samples and it cannot invent what is not
# in the file -- most importantly the grouping, without which there is
# nothing to run a differential analysis between, and the donor, without
# which two layers cannot be paired.
#
# A template is cheaper than documenting the rules: the shape is the
# documentation, and a file that already parses is a much better starting
# point than a description of one.
#
# The two templates deliberately use *different* sample ids and the *same*
# donors. That is the whole integration story in one glance: the same six
# people, measured twice, and the column that says so.

TEMPLATE_DONORS <- sprintf("D%02d", 1:6)
TEMPLATE_CONDITION <- c("G1", "G1", "G1", "G2", "G2", "G2")

#' Sheets making up an import template
#'
#' @param omics_type `"proteomics"` or `"rnaseq"`.
#' @return Named list of data frames, in the order they should be written.
#' @keywords internal
#' @noRd
import_template_sheets <- function(omics_type = c("proteomics", "rnaseq")) {
  omics_type <- match.arg(omics_type)
  set.seed(1L)

  n_feat <- 8L
  sample_ids <- if (omics_type == "proteomics") {
    sprintf("P%02d", seq_along(TEMPLATE_DONORS))
  } else {
    sprintf("R%02d", seq_along(TEMPLATE_DONORS))
  }

  # Real ids, not plausible-looking ones: the RNA-seq template exists
  # partly to show that symbols arrive by mapping, and a made-up ENSG
  # would map to nothing and teach the opposite.
  feature_ids <- if (omics_type == "proteomics") {
    c("P01308", "P02768", "P69905", "P68871",
      "P01009", "P02787", "P00738", "P01023")
  } else {
    c("ENSG00000254647", "ENSG00000163631", "ENSG00000206172",
      "ENSG00000244734", "ENSG00000197249", "ENSG00000091513",
      "ENSG00000257017", "ENSG00000175899")
  }
  symbols <- c("INS", "ALB", "HBA1", "HBB",
               "SERPINA1", "TF", "HP", "A2M")

  values <- if (omics_type == "proteomics") {
    # Intensities: continuous, wide dynamic range.
    round(matrix(stats::rlnorm(n_feat * length(sample_ids), 12, 1.2),
                 nrow = n_feat), 1)
  } else {
    # Counts: integers, and some genes genuinely zero everywhere.
    m <- matrix(stats::rnbinom(n_feat * length(sample_ids), mu = 200, size = 2),
                nrow = n_feat)
    m[n_feat, ] <- 0L
    m
  }

  expression <- data.frame(feature_id = feature_ids, values,
                           check.names = FALSE, stringsAsFactors = FALSE)
  names(expression)[-1L] <- sample_ids

  sample_info <- data.frame(
    sample_id = sample_ids,
    donor     = TEMPLATE_DONORS,
    condition = TEMPLATE_CONDITION,
    age       = c(24, 31, 28, 55, 61, 58),
    sex       = c("F", "M", "F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  out <- list(
    README      = import_template_readme(omics_type, sample_ids),
    expression  = expression,
    sample_info = sample_info
  )

  # Two data sheets for RNA-seq, three for proteomics -- each mirrors what
  # that kind of vendor actually ships, rather than forcing one shape on
  # both.
  #
  # RNA-seq needs no feature sheet because the symbol is derived from the
  # Ensembl id against a current HGNC table, which is better than any
  # symbol column a vendor froze on the day they ran the pipeline.
  # Proteomics keeps one: a UniProt accession has no such mapping here,
  # so the symbol has to come from the file.
  if (omics_type == "proteomics") {
    out$feature_info <- data.frame(
      feature_id          = feature_ids,
      gene_symbol         = symbols,
      protein_description = paste(symbols, "protein"),
      stringsAsFactors    = FALSE
    )
  }
  out
}

# Written as a one-column sheet so it travels inside the workbook. The
# classifier gives it role "unknown" and ignores it, which is the right
# outcome and visible in the Sheet roles list.
import_template_readme <- function(omics_type, sample_ids) {
  values <- if (omics_type == "proteomics") {
    "Intensities. Any scale; the import asks you which."
  } else {
    "Raw counts, integers. Not FPKM or TPM -- DESeq2 needs counts."
  }

  data.frame(
    `How to fill this in` = c(
      sprintf("Template for a %s layer. Replace the example rows.", omics_type),
      "",
      "SHEET 'expression'",
      "  Column 1 is the feature id. Every other column is one sample.",
      paste0("  Values: ", values),
      "",
      "SHEET 'sample_info'",
      "  One row per sample.",
      "  sample_id  must match the expression column headings exactly.",
      "  condition  the grouping compared in Differential. Without it",
      "             there is nothing to compare and the view stays empty.",
      "  donor      the person. This is what pairs two layers in",
      "             Integration -- the proteomics and RNA-seq templates",
      "             use different sample_ids and the SAME donors, which",
      "             is what a real study looks like: one person, two",
      "             tissues, two sample names.",
      "             Optional. Filling it in here means Integration is",
      "             already paired up; leaving it out means pairing the",
      "             samples in the Integration view instead. Nothing has",
      "             to be re-imported or re-run either way -- donor is",
      "             read by Integration and by nothing else.",
      "  Extra columns are kept and can be used for grouping or colour.",
      "",
      symbol_section(omics_type),
      "",
      "The sheet NAMES do not matter -- the importer classifies by",
      "content and lets you correct it. The COLUMN CONTENTS do.",
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Where the gene symbol comes from is the one thing that differs between
# the two templates, and it is the thing that decides whether enrichment
# returns anything -- so it gets said explicitly rather than left to be
# inferred from which sheets are present.
symbol_section <- function(omics_type) {
  if (identical(omics_type, "rnaseq")) {
    c("GENE SYMBOLS -- nothing to fill in",
      "  Two sheets is all an RNA-seq layer needs. Symbols are looked up",
      "  from the Ensembl id against a current HGNC table, which beats a",
      "  gene_name column frozen on the day the pipeline ran: symbols get",
      "  renamed, merged and retired, ids do not.",
      "  The import report says how many matched. Ids with no HGNC symbol",
      "  -- unnamed lncRNAs, pseudogenes -- are still tested in",
      "  Differential and are left out of enrichment, which no pathway",
      "  database would have matched anyway.")
  } else {
    c("SHEET 'feature_info'",
      "  feature_id   must match column 1 of 'expression'.",
      "  gene_symbol  what enrichment matches on. A UniProt accession is",
      "               not a gene name, and there is no id mapping for",
      "               proteomics here -- so without this column every",
      "               database returns nothing.")
  }
}

#' Write an import template workbook
#'
#' @param path Destination `.xlsx`.
#' @param omics_type `"proteomics"` or `"rnaseq"`.
#' @return `path`, invisibly.
#' @keywords internal
#' @noRd
write_import_template <- function(path, omics_type = "proteomics") {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required to write the template workbook.",
         call. = FALSE)
  }
  sheets <- import_template_sheets(omics_type)
  openxlsx::write.xlsx(sheets, file = path, overwrite = TRUE)
  invisible(path)
}
