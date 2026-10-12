# Gzip-compressed text: counts.csv.gz, quant.sf.gz, a GEO supplementary
# table. Unpacked to a temporary file first, so every reader downstream
# -- the encoding check, fread, read.table, the quantification-file
# check -- sees the plain text it already handles.
#
# The unpacking is where a small upload turns into a large one. A gzip
# stream records its unpacked size in its last four bytes, and
# guard_archive() refuses a file whose record is over the limit before
# anything is read; but the record is only the size modulo 4 GB, and a
# file made of several gzip members records only the last. So the bytes
# are also counted as they come out, and unpacking stops at the same
# limit -- `options(omicsCore.max_unpacked_mb =)`, as for workbooks.

GZIP_MAGIC <- as.raw(c(0x1f, 0x8b))

# Whether the file starts like a gzip stream. The bytes decide, not the
# name: a browser upload arrives as "0.gz" or, renamed, as anything.
is_gzip_file <- function(path) {
  con <- file(path, "rb")
  magic <- tryCatch(readBin(con, "raw", 2L), finally = close(con))
  length(magic) == 2L && identical(magic, GZIP_MAGIC)
}

#' The plain text of a possibly gzip-compressed file
#'
#' @param path Path to a text file, compressed or not.
#' @param max_mb Most megabytes the unpacked text may take.
#' @return `list(path, converted)`. `path` is `path` itself when it was
#'   not compressed, and otherwise a temporary file the caller removes
#'   (`converted = TRUE`).
#' @keywords internal
#' @noRd
gunzip_text_file <- function(path,
                             max_mb = getOption("omicsCore.max_unpacked_mb", 2048)) {
  unchanged <- list(path = path, converted = FALSE)
  if (!is_gzip_file(path)) return(unchanged)
  limit <- max_mb * 1024^2
  out <- tempfile("omics-gunzip-", fileext = ".txt")
  src <- gzfile(path, "rb")
  dst <- file(out, "wb")
  written <- 0
  status <- tryCatch({
    repeat {
      chunk <- readBin(src, "raw", n = 8L * 1024L^2)
      if (!length(chunk)) break
      written <- written + length(chunk)
      if (written > limit) {
        stop(structure(class = c("omics_gunzip_limit", "error", "condition"),
                       list(message = "", call = NULL)))
      }
      writeBin(chunk, dst)
    }
    "ok"
  },
  omics_gunzip_limit = function(e) "too_big",
  # A truncated or damaged stream: zlib's complaint is about zlib.
  error = function(e) "damaged",
  warning = function(w) "damaged",
  finally = {
    close(src)
    close(dst)
  })
  if (!identical(status, "ok")) {
    unlink(out)
    if (identical(status, "too_big")) {
      stop(sprintf(paste(
        "This compressed file unpacks to more than the %s this server reads.",
        "Upload a smaller file, or raise options(omicsCore.max_unpacked_mb)."),
        format_mb(limit)), call. = FALSE)
    }
    stop("The compressed file could not be unpacked; it may be damaged or ",
         "incomplete. Unpack it on your computer and upload the file inside.",
         call. = FALSE)
  }
  list(path = out, converted = TRUE)
}

# Extensions of the text a .gz may hold, as read_omics() names them.
GZ_TEXT_EXTENSIONS <- c("csv", "tsv", "txt", "tab", "sf", "results")
