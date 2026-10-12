# Text files that are not UTF-8.
#
# A CSV saved by Excel on Windows is written in the machine's "ANSI" code
# page, not in UTF-8: Windows-1252 in Western Europe and the Americas,
# GBK (the two-byte part of GB18030) on Chinese Windows. Read as UTF-8,
# the first "µg/ml" in a unit column or a Chinese group name is an
# invalid string, and the import stopped with "input string 2 is invalid
# UTF-8" -- a message about R, not about the file. The bytes are looked
# at first and the file is decoded to UTF-8 before anything parses it,
# so every reader downstream sees the UTF-8 it already expects.
#
# The order of the guesses:
#   * a UTF-16 byte-order mark: Excel's "Unicode text" export;
#   * plain ASCII, or valid UTF-8 (with or without a BOM): nothing to do;
#   * valid GB18030 whose high bytes come in pairs that both lie in
#     0xA1-0xFE, the shape of GB2312's common characters. Validity alone
#     does not separate it from Windows-1252: "µg" (B5 67) and "°C"
#     (B0 43) are valid GBK pairs too, decoding to rare ideographs. In
#     Western text an accented letter stands alone between ASCII
#     letters; in Chinese text the high bytes come in long even runs;
#   * otherwise Windows-1252, or Latin-1 when the file holds one of the
#     five bytes Windows-1252 leaves undefined (Latin-1 defines all 256,
#     so it always decodes).

# Plain-English names for the report.
TEXT_ENCODING_LABELS <- c(
  `UTF-16` = "UTF-16",
  GB18030 = "Chinese GB18030/GBK",
  CP1252 = "Windows-1252",
  latin1 = "Latin-1"
)

# Bytes Windows-1252 does not define.
CP1252_UNDEFINED <- as.raw(c(0x81, 0x8d, 0x8f, 0x90, 0x9d))

# Calls `fun(buf)` on successive slices of the file, each ending at a
# newline, so a multi-byte character is never cut in two (0x0A is never
# part of one in UTF-8, GBK or the single-byte code pages). The slices
# keep a large file from being held in memory whole twice over.
for_each_text_chunk <- function(path, fun, chunk_bytes = 8 * 1024^2) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  carry <- raw(0)
  repeat {
    chunk <- readBin(con, "raw", n = chunk_bytes)
    done <- length(chunk) == 0L
    buf <- c(carry, chunk)
    carry <- raw(0)
    if (!done && length(buf)) {
      # The last line end, looked for near the end first: scanning the
      # whole slice for every newline costs more than the rest together.
      n <- length(buf)
      from <- max(1L, n - 65535L)
      nl <- which(buf[from:n] == as.raw(0x0a))
      cut <- if (length(nl)) from - 1L + nl[[length(nl)]] else {
        nl <- which(buf == as.raw(0x0a))
        if (length(nl)) nl[[length(nl)]] else 0L
      }
      if (cut < n) carry <- buf[(cut + 1L):n]
      # readBin() on the vector: a head without the copy subsetting makes.
      buf <- readBin(buf, "raw", cut)
    }
    # `fun` returns FALSE once it has seen enough.
    if (length(buf) && identical(fun(buf), FALSE)) break
    if (done) break
  }
  invisible(NULL)
}

# Whether any byte is outside ASCII.
file_has_high_byte <- function(path, chunk_bytes = 8 * 1024^2) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  repeat {
    chunk <- readBin(con, "raw", n = chunk_bytes)
    if (!length(chunk)) return(FALSE)
    if (any(chunk >= as.raw(0x80))) return(TRUE)
  }
}

#' Which text encoding a file is in
#'
#' @param path Path to a delimited text file.
#' @return One of `"UTF-8"` (which covers plain ASCII), `"UTF-16"`,
#'   `"GB18030"`, `"CP1252"` or `"latin1"`.
#' @keywords internal
#' @noRd
detect_text_encoding <- function(path, chunk_bytes = 8 * 1024^2) {
  con <- file(path, "rb")
  head <- tryCatch(readBin(con, "raw", 2L), finally = close(con))
  if (length(head) == 2L &&
      (identical(head, as.raw(c(0xff, 0xfe))) || identical(head, as.raw(c(0xfe, 0xff))))) {
    return("UTF-16")
  }
  # Most files are plain ASCII, and this is all they cost: one look at
  # the bytes, without cutting them at line ends.
  if (!file_has_high_byte(path, chunk_bytes)) return("UTF-8")

  # Then UTF-8, which is the rest of them, stopping at the first slice
  # that is not. A NUL byte is not text in any of these encodings; the
  # readers report that file as they always have.
  utf8_ok <- TRUE
  binary <- FALSE
  for_each_text_chunk(path, chunk_bytes = chunk_bytes, function(buf) {
    # rawToChar() refuses an embedded NUL, which is the check.
    s <- tryCatch(rawToChar(buf), error = function(e) NULL)
    if (is.null(s)) {
      binary <<- TRUE
      return(FALSE)
    }
    if (!validUTF8(s)) {
      utf8_ok <<- FALSE
      return(FALSE)
    }
    TRUE
  })
  if (binary || utf8_ok) return("UTF-8")

  # Not UTF-8: GB18030 or a Western code page.
  gb_ok <- TRUE
  cp_ok <- TRUE
  n_high <- 0
  n_gb_shaped <- 0
  for_each_text_chunk(path, chunk_bytes = chunk_bytes, function(buf) {
    hi <- which(buf >= as.raw(0x80))
    if (!length(hi)) return(TRUE)
    if (gb_ok && is.na(iconv(rawToChar(buf), "GB18030", "UTF-8"))) gb_ok <<- FALSE
    if (cp_ok && any(buf[hi] %in% CP1252_UNDEFINED)) cp_ok <<- FALSE
    # Runs of consecutive high bytes; a run is GB2312-shaped when it is
    # of even length and none of its bytes is below 0xA1.
    run <- cumsum(c(TRUE, diff(hi) != 1L))
    run_len <- tabulate(run)
    low <- buf[hi] < as.raw(0xa1)
    has_low <- tabulate(run[low], nbins = length(run_len)) > 0L
    shaped <- run_len %% 2L == 0L & !has_low
    n_high <<- n_high + length(hi)
    n_gb_shaped <<- n_gb_shaped + sum(run_len[shaped])
    TRUE
  })
  if (gb_ok && n_gb_shaped >= 0.8 * n_high) return("GB18030")
  if (cp_ok) "CP1252" else "latin1"
}

#' A UTF-8 copy of a text file, when it is not UTF-8 already
#'
#' @param path Path to a delimited text file.
#' @return `list(path, encoding, converted)`. `path` is `path` itself
#'   when nothing had to change, and otherwise a temporary UTF-8 file the
#'   caller removes. `encoding` is what the file was read as.
#' @keywords internal
#' @noRd
utf8_text_file <- function(path) {
  enc <- tryCatch(detect_text_encoding(path), error = function(e) "UTF-8")
  unchanged <- list(path = path, encoding = "UTF-8", converted = FALSE)
  if (identical(enc, "UTF-8")) return(unchanged)
  out <- tempfile("omics-utf8-", fileext = ".txt")
  ok <- tryCatch({
    if (identical(enc, "UTF-16")) {
      bytes <- readBin(path, "raw", n = file.size(path))
      # iconv() on raw input hands back the input unchanged when it
      # cannot convert it, so success is judged on the result: UTF-8
      # text holds no NUL bytes, and UTF-16 of ASCII is half NULs.
      conv <- iconv(list(bytes), "UTF-16", "UTF-8", toRaw = TRUE)[[1L]]
      if (is.null(conv) || any(conv == as.raw(0L))) stop("not UTF-16")
      writeBin(conv, out)
    } else {
      con <- file(out, "wb")
      on.exit(close(con), add = TRUE)
      for_each_text_chunk(path, function(buf) {
        s <- iconv(rawToChar(buf), enc, "UTF-8")
        # Windows-1252 was chosen because it defines every byte present,
        # and Latin-1 defines all of them; a failure here is a file that
        # changed under us.
        if (is.na(s)) s <- iconv(rawToChar(buf), "latin1", "UTF-8")
        writeBin(charToRaw(s), con)
      })
    }
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    unlink(out)
    return(unchanged)
  }
  list(path = out, encoding = enc, converted = TRUE)
}

# The note the import page shows. NULL for UTF-8, which is the
# expectation and needs no comment.
encoding_note <- function(encoding, what = "The file") {
  if (is.null(encoding) || identical(encoding, "UTF-8")) return(NULL)
  sprintf("%s was read as %s (not UTF-8).", what,
          unname(TEXT_ENCODING_LABELS[encoding]) %|NA|% encoding)
}

`%|NA|%` <- function(a, b) if (length(a) == 0L || is.na(a)) b else a
