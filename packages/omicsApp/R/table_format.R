# Numbers in the on-screen result tables.
#
# A p-value went to the browser as a number, and the browser printed
# it in full: signif(1.8e-248, 3) is not exactly representable, and the
# table said "1.799999999999999e-248" -- 22 characters with no place to
# break, which pushed the last column out of the card. The tables now
# send the p-value as text, written the way the rest of the app writes
# it, and sort that column by a hidden copy of the number: as text,
# "9.1e-05" would sort above "0.012".
#
# The files a table downloads are written from the result itself, at
# full precision; none of this touches them.

# 0.0123, 0.5, 4.42e-16: three significant digits, scientific below
# 0.001. Vectorised; a missing value is a dash.
format_p_value <- function(p, digits = 3L) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep("\u2013", length(p))
  ok <- !is.na(p)
  # A p-value that underflowed to zero is smaller than any written one.
  zero <- ok & p <= 0
  out[zero] <- "0"
  # Rounded first, so 0.0009999 reads 0.001 like 0.001 itself.
  small <- ok & !zero & signif(p, digits) < 1e-3
  out[small] <- formatC(p[small], format = "e", digits = digits - 1L)
  # One at a time: format() pads a vector to its widest member.
  plain <- ok & !zero & !small
  out[plain] <- vapply(p[plain], function(v) format(signif(v, digits), scientific = FALSE),
                       character(1))
  out
}

#' A results table's numbers as text that still sorts as numbers
#'
#' Replaces each of `p_cols` with format_p_value() text, and each column
#' named in `sort_keys` keeps its text but sorts by the number given.
#' Both are right-aligned.
#' The numbers travel as hidden columns at the end of the table, so a
#' selected row's index and the visible column positions are unchanged.
#'
#' @param out The data frame as it will be shown.
#' @param p_cols Names of its p-value columns (numeric).
#' @param sort_keys Named list: shown column name -> numeric vector to
#'   sort it by.
#' @return list(data = the frame to hand to DT::datatable(),
#'   column_defs = the columnDefs entries that hide and use the keys).
#' @keywords internal
#' @noRd
dt_sortable_text <- function(out, p_cols = character(0), sort_keys = list()) {
  keys <- list()
  for (col in intersect(p_cols, names(out))) {
    keys[[col]] <- suppressWarnings(as.numeric(out[[col]]))
    out[[col]] <- format_p_value(out[[col]])
  }
  for (col in intersect(names(sort_keys), names(out))) keys[[col]] <- sort_keys[[col]]
  defs <- list()
  if (length(keys)) {
    n <- ncol(out)
    for (i in seq_along(keys)) {
      out[[paste0(".sort_", i)]] <- as.numeric(keys[[i]])
      # DataTables counts columns from 0. A number, so on the right.
      defs[[length(defs) + 1L]] <- list(targets = match(names(keys)[[i]], names(out)) - 1L,
                                        orderData = n + i - 1L, className = "dt-right")
    }
    defs[[length(defs) + 1L]] <- list(targets = n + seq_along(keys) - 1L,
                                      visible = FALSE, searchable = FALSE)
  }
  list(data = out, column_defs = defs)
}
