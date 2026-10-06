# Atomic save / load of an omics_project. We use qs2 because it round-trips
# arbitrary R objects (DESeqDataSet, MArrayLM, ComplexHeatmap, ggplot, etc.)
# losslessly and is faster than saveRDS on the kinds of objects analysis
# bundles tend to hold. qs2 is Suggests-gated so a base install of
# omicsCore stays slim.
#
# File format: a qs2 archive whose root is a list with three fields:
#   * `schema_version`  -- semantic-version-ish string. Kept at "1.x" for
#                          as long as older readers can still open the
#                          file: they refuse anything whose major part
#                          is higher, so this is the switch that keeps a
#                          rolled-back release away from a file it would
#                          misread.
#   * `format_version`  -- integer, the version the migrations below are
#                          keyed by. Files written before it existed
#                          carry none and are format 1.
#   * `payload`         -- the `omics_project`
#
# Optionally followed by a signature trailer (see "Signing" below), which
# qs2 ignores: a release that predates signing still reads a signed file,
# which is what makes rolling back an image safe.

OMP_SCHEMA_VERSION <- "1.0.0"

# Bump when the payload's layout changes, and register a migration from
# the previous version (register_omp_migration()) in the same commit.
OMP_FORMAT_VERSION <- 1L

# A function rather than a bare constant so a test can pretend the
# current version is newer and exercise a migration end to end.
current_omp_format_version <- function() OMP_FORMAT_VERSION

ensure_qs2 <- function() {
  if (!is_installed("qs2")) {
    stop(
      "Package 'qs2' is required for save_project() / load_project(). ",
      "Install with: install_optional('persistence').",
      call. = FALSE
    )
  }
}

# ---- migrations ------------------------------------------------------------
#
# One function per step, keyed by the version it upgrades *from*: the
# function registered under "1" takes a format-1 project and returns a
# format-2 one. load_project() chains them, so a file several versions
# old walks through each step in turn and no migration ever has to know
# about more than one change.

.omp_migrations <- new.env(parent = emptyenv())

#' Register a migration between two `.omp` format versions
#'
#' @param from Integer format version the migration upgrades from.
#' @param fn Function taking the project as stored in format `from` and
#'   returning it in format `from + 1`.
#' @return Invisibly, the migration previously registered for `from`
#'   (or `NULL`).
#' @keywords internal
#' @noRd
register_omp_migration <- function(from, fn) {
  if (!is.numeric(from) || length(from) != 1L || is.na(from) ||
      from < 1 || from != round(from)) {
    stop("`from` must be a single whole number of at least 1.")
  }
  if (!is.function(fn)) stop("`fn` must be a function.")
  key <- as.character(as.integer(from))
  old <- if (exists(key, envir = .omp_migrations, inherits = FALSE)) {
    get(key, envir = .omp_migrations)
  }
  assign(key, fn, envir = .omp_migrations)
  invisible(old)
}

# The project in `envelope`, upgraded to `target`.
migrate_omp_payload <- function(envelope, target = current_omp_format_version()) {
  from <- envelope$format_version
  if (is.null(from)) from <- 1L
  if (!is.numeric(from) || length(from) != 1L || is.na(from) ||
      from < 1 || from != round(from)) {
    stop("Invalid format_version in the project file.")
  }
  from <- as.integer(from)
  if (from > target) {
    stop("This project was saved by a newer version of the app (file format ",
         from, "); this version reads up to format ", target,
         ". Update the app to open it.")
  }
  project <- envelope$payload
  while (from < target) {
    key <- as.character(from)
    if (!exists(key, envir = .omp_migrations, inherits = FALSE)) {
      stop("No upgrade is available for project file format ", from, ".")
    }
    project <- get(key, envir = .omp_migrations)(project)
    from <- from + 1L
  }
  project
}

# ---- signing ---------------------------------------------------------------
#
# A signed file is the qs2 archive followed by a 40-byte trailer: the
# HMAC-SHA256 of the archive under the caller's key, then the 8-byte
# magic below. The magic goes last so a reader can tell a signed file
# from an unsigned one by its final bytes alone, before deserialising
# anything.
#
# The point is that deserialising is the risky step. A qs2 file can hold
# any R object, and a reader that checks a file's shape only after
# reading it has already let the file decide what to build in memory. A
# file that verifies under a key only this deployment holds was written
# by this deployment, so it is read; anything else is refused before qs2
# sees a byte of it.

OMP_SIGNATURE_MAGIC <- as.raw(c(0x4f, 0x4d, 0x50, 0x53, 0x49, 0x47, 0x30, 0x31)) # "OMPSIG01"
OMP_MAC_BYTES <- 32L
OMP_TRAILER_BYTES <- OMP_MAC_BYTES + length(OMP_SIGNATURE_MAGIC)

ensure_openssl <- function() {
  if (!is_installed("openssl")) {
    stop("Package 'openssl' is required to sign or verify project files. ",
         "Install with: install.packages(\"openssl\").", call. = FALSE)
  }
}

# A key as raw bytes, or NULL. A string is taken as its UTF-8 bytes, so
# the same text in an environment variable and in a key file is the same
# key.
as_signing_key <- function(key, arg = "signing_key") {
  if (is.null(key)) return(NULL)
  if (is.raw(key) && length(key) > 0L) return(key)
  if (is.character(key) && length(key) == 1L && !is.na(key) && nzchar(key)) {
    return(charToRaw(enc2utf8(key)))
  }
  stop("`", arg, "` must be NULL, or a non-empty string or raw vector.")
}

omp_mac <- function(body, key) {
  ensure_openssl()
  as.raw(openssl::sha256(body, key = key))
}

# Compares every byte whatever the first difference, so the time taken
# says nothing about how much of a forged signature was right.
mac_equal <- function(a, b) {
  length(a) == length(b) && sum(as.integer(xor(a, b))) == 0L
}

# list(body, mac): `mac` is NULL when the bytes carry no signature trailer.
split_omp_signature <- function(bytes) {
  n <- length(bytes)
  if (n <= OMP_TRAILER_BYTES ||
      !identical(bytes[(n - length(OMP_SIGNATURE_MAGIC) + 1L):n], OMP_SIGNATURE_MAGIC)) {
    return(list(body = bytes, mac = NULL))
  }
  body_n <- n - OMP_TRAILER_BYTES
  list(body = bytes[seq_len(body_n)],
       mac = bytes[(body_n + 1L):(body_n + OMP_MAC_BYTES)])
}

read_file_bytes <- function(path) {
  size <- file.size(path)
  if (is.na(size) || size == 0) return(raw(0))
  readBin(path, "raw", n = size)
}

omp_condition <- function(subclass, message) {
  structure(class = c(subclass, "omp_untrusted_error", "error", "condition"),
            list(message = message, call = NULL))
}

# The file's archive bytes, verified under `key`. Errors (classed
# `omp_untrusted_error`) when the file is unsigned or the signature does
# not match.
verified_omp_body <- function(bytes, key) {
  parts <- split_omp_signature(bytes)
  if (is.null(parts$mac)) {
    stop(omp_condition("omp_unsigned_error",
      "This file was not saved by this app, so it was not opened."))
  }
  if (!mac_equal(omp_mac(parts$body, key), parts$mac)) {
    # Also what a file signed under some other key gets: to the user the
    # two are the same -- this app did not write these bytes.
    stop(omp_condition("omp_signature_error",
      "This file has been changed since this app saved it, so it was not opened."))
  }
  parts$body
}

# Writes `body` plus a signature under `key` to `path`, atomically.
write_signed_omp <- function(body, key, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."),
                  tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
  con <- file(tmp, open = "wb")
  ok <- FALSE
  tryCatch({
    writeBin(body, con)
    writeBin(c(omp_mac(body, key), OMP_SIGNATURE_MAGIC), con)
    ok <- TRUE
  }, finally = close(con))
  if (!isTRUE(ok) || !isTRUE(file.rename(tmp, path))) {
    stop("Failed to write the signed file '", path, "'.")
  }
  invisible(path)
}

# ---- what an untrusted file may contain -------------------------------------
#
# A file the app did not sign is checked after it is read and before
# anything uses it. A project is data: numbers, strings, lists, data
# frames, and the result objects of the analysis packages. A function or
# an environment in one is code that would run with the user's
# permissions the moment something called or printed it, and an
# external pointer is a C-level address no file has any business
# carrying -- so all three are refused wherever they sit, as are classes
# the app never writes.

OMP_SAFE_TYPES <- c("NULL", "logical", "integer", "double", "complex",
                    "character", "raw", "list", "S4")

OMP_SAFE_CLASSES <- c(
  # omicsCore
  "omics_project", "omics_input", "analysis_bundle", "ImportReport",
  # base R data
  "data.frame", "tbl_df", "tbl", "factor", "ordered", "Date", "POSIXct",
  "POSIXt", "difftime", "table", "AsIs", "noquote", "numeric_version",
  "package_version", "matrix", "array",
  # stats / limma / edgeR result objects kept in bundles
  "qr", "MArrayLM", "TestResults", "DGEList", "DGEGLM", "DGELRT", "DGEExact",
  # clusterProfiler / DOSE enrichment objects
  "enrichResult", "gseaResult", "compareClusterResult"
)

OMP_MAX_DEPTH <- 64L

describe_type <- function(type) {
  switch(type,
         closure = , builtin = , special = "a function",
         environment = "an environment",
         externalptr = "an external pointer",
         promise = "a promise (code waiting to run)",
         language = , symbol = , expression = "R code",
         paste0("an object of type '", type, "'"))
}

#' Check that an object read from an untrusted file holds only data
#'
#' Walks every element and every attribute. Errors (class
#' `omp_unsafe_error`) naming the first thing that is not allowed.
#'
#' @param x The deserialised object.
#' @return Invisibly `TRUE`.
#' @keywords internal
#' @noRd
check_project_structure <- function(x) {
  refuse <- function(what, where) {
    stop(omp_condition("omp_unsafe_error", sprintf(
      "This file contains %s (at %s), which a project never holds, so it was not opened.",
      what, where)))
  }
  walk <- function(obj, where, depth) {
    if (depth > OMP_MAX_DEPTH) refuse("data nested too deeply", where)
    type <- typeof(obj)
    if (!type %in% OMP_SAFE_TYPES) refuse(describe_type(type), where)
    cls <- attr(obj, "class", exact = TRUE)
    if (!is.null(cls)) {
      bad <- setdiff(as.character(cls), OMP_SAFE_CLASSES)
      if (length(bad)) refuse(sprintf("an object of class '%s'", bad[[1L]]), where)
    }
    # Every attribute, structural ones included: S4 slots are attributes,
    # and `names` or `dimnames` can hold anything a file says they do.
    attrs <- attributes(obj)
    for (nm in names(attrs)) {
      walk(attrs[[nm]], paste0(where, "@", nm), depth + 1L)
    }
    if (type == "list") {
      nms <- names(obj)
      for (i in seq_along(obj)) {
        label <- if (!is.null(nms) && !is.na(nms[[i]]) && nzchar(nms[[i]])) {
          paste0("$", nms[[i]])
        } else {
          paste0("[[", i, "]]")
        }
        walk(obj[[i]], paste0(where, label), depth + 1L)
      }
    }
    invisible(TRUE)
  }
  walk(x, "the top level", 0L)
}

# From a deserialised envelope to a project: shape, version, migrations.
project_from_envelope <- function(envelope, path, untrusted = FALSE) {
  if (isTRUE(untrusted)) check_project_structure(envelope)
  if (!is.list(envelope) ||
      !all(c("schema_version", "payload") %in% names(envelope))) {
    stop("`", path, "` does not appear to be an omicsCore project file.")
  }
  schema <- envelope$schema_version
  if (!is.character(schema) || length(schema) != 1L) {
    stop("Invalid schema_version in '", path, "'.")
  }
  if (!schema_is_supported(schema)) {
    stop("Unsupported project schema version: '", schema,
         "'. This omicsCore supports up to '", OMP_SCHEMA_VERSION, "'.")
  }
  project <- migrate_omp_payload(envelope)
  if (!is_omics_project(project)) {
    stop("Loaded payload is not an `omics_project`.")
  }
  project
}

#' Save an omics_project to disk
#'
#' Writes the full project (every experiment, sample_link, and any
#' attached `analysis_bundle`s) to a single `.omp` file using `qs2`. The
#' write is atomic: the payload first goes to `path.tmp`, then renames
#' over `path`, so an interrupted save will not corrupt an existing file.
#'
#' With `signing_key`, the file ends in an HMAC-SHA256 signature under
#' that key, which [load_project()] checks before reading anything when
#' given the same key. Readers that do not ask for a signature ignore it.
#'
#' Requires the `qs2` package; install it with `install_optional("persistence")`.
#' Signing also requires `openssl`.
#'
#' @param project An [`omics_project`][is_omics_project()].
#' @param path Output path. The file extension is up to the caller, but
#'   the convention is `.omp` for projects produced by `omicsCore`.
#' @param overwrite If `FALSE` (default) and `path` already exists,
#'   `save_project()` raises an error. Set to `TRUE` to replace.
#' @param signing_key `NULL` (default) for an unsigned file, or the key to
#'   sign it with: a non-empty string (taken as its UTF-8 bytes) or a raw
#'   vector.
#'
#' @return Invisibly returns `path`.
#' @export
#' @family persistence
#' @examples
#' \dontrun{
#'   p <- omics_project("demo", experiments = list(proteo = my_input))
#'   save_project(p, "demo.omp")
#'   q <- load_project("demo.omp")
#' }
save_project <- function(project, path, overwrite = FALSE, signing_key = NULL) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty single string.")
  }
  key <- as_signing_key(signing_key)
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Path already exists: ", path,
         " (pass `overwrite = TRUE` to replace).")
  }
  ensure_qs2()
  if (!is.null(key)) ensure_openssl()

  payload <- list(
    schema_version = OMP_SCHEMA_VERSION,
    format_version = current_omp_format_version(),
    payload = project
  )

  # Unique per call, not paste0(path, ".tmp").
  #
  # The rename is what makes this atomic for a reader, and a shared
  # temp name takes that back: two writers -- the same account open on
  # two devices, both autosaving -- open the same file, interleave, and
  # one renames the result over the real one. The on.exit cleanup would
  # also delete the other writer's file mid-write.
  #
  # Same directory as the target, because rename() is only atomic within
  # a filesystem and tempdir() is usually a different one.
  tmp <- tempfile(pattern = paste0(basename(path), "."),
                  tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
  qs2::qs_save(payload, file = tmp)
  if (!is.null(key)) {
    # Appended to the temporary file, so the signature lands in the same
    # rename as the archive: no reader ever sees one without the other.
    mac <- omp_mac(read_file_bytes(tmp), key)
    con <- file(tmp, open = "ab")
    tryCatch(writeBin(c(mac, OMP_SIGNATURE_MAGIC), con), finally = close(con))
  }
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    stop("Failed to rename '", tmp, "' to '", path, "'.")
  }
  invisible(path)
}

#' Load an omics_project from disk
#'
#' Reads an `.omp` archive written by [save_project()]. Validates the
#' schema version and that the payload is an `omics_project`, and
#' upgrades a file written in an older format through the registered
#' migrations.
#'
#' With `signing_key`, the file must carry a valid signature under that
#' key; an unsigned or altered file is refused before it is deserialised.
#' With `untrusted = TRUE`, the object read is checked to hold only data
#' -- no functions, environments, external pointers or unexpected classes
#' -- before it is returned. Use it for a file that came from anywhere
#' other than your own store.
#'
#' Requires the `qs2` package; install it with `install_optional("persistence")`.
#' Verifying a signature also requires `openssl`.
#'
#' @param path Path to the `.omp` file.
#' @param signing_key `NULL` (default) to read without checking a
#'   signature, or the key the file must be signed with (see
#'   [save_project()]).
#' @param untrusted If `TRUE`, check the structure of what was read
#'   before returning it. Default `FALSE`.
#'
#' @return The deserialized `omics_project`.
#' @export
#' @family persistence
load_project <- function(path, signing_key = NULL, untrusted = FALSE) {
  if (!is.character(path) || length(path) != 1L) {
    stop("`path` must be a single string.")
  }
  key <- as_signing_key(signing_key)
  if (!is.logical(untrusted) || length(untrusted) != 1L || is.na(untrusted)) {
    stop("`untrusted` must be TRUE or FALSE.")
  }
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  ensure_qs2()
  envelope <- if (is.null(key)) {
    qs2::qs_read(path)
  } else {
    # Verified as a separate step, before qs2 is handed anything.
    body <- verified_omp_body(read_file_bytes(path), key)
    qs2::qs_deserialize(body)
  }
  project_from_envelope(envelope, path, untrusted = untrusted)
}

#' Sign a project file in place
#'
#' Used by a project store to adopt files it did not write itself, and to
#' move its files to a new key. The file is read once; what is signed is
#' exactly what was checked.
#'
#' * A file already signed under `signing_key` is left alone.
#' * A file signed under `verify_key` (a previous key) is re-signed under
#'   `signing_key`.
#' * An unsigned file is adopted only when `adopt_unsigned = TRUE`, and
#'   only if it reads as a project holding nothing but data (see
#'   [load_project()]'s `untrusted`).
#' * Anything else -- a file whose signature matches neither key, or an
#'   unsigned file that is not adopted -- is left unchanged.
#'
#' The archive bytes are never rewritten, so a signed file reads exactly
#' as it did before.
#'
#' @param path Path to an `.omp` file.
#' @param signing_key The key to sign with (see [save_project()]).
#' @param verify_key Optional previous key whose signatures are accepted
#'   and replaced.
#' @param adopt_unsigned Whether to sign an unsigned file that passes the
#'   structure check.
#'
#' @return A list with `status` -- one of `"valid"`, `"signed"`,
#'   `"resigned"` or `"refused"` -- and `message`.
#' @export
#' @family persistence
sign_project_file <- function(path, signing_key, verify_key = NULL,
                              adopt_unsigned = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a non-empty single string.")
  }
  key <- as_signing_key(signing_key)
  if (is.null(key)) stop("`signing_key` must be a non-empty string or raw vector.")
  old <- as_signing_key(verify_key, "verify_key")
  if (!is.logical(adopt_unsigned) || length(adopt_unsigned) != 1L || is.na(adopt_unsigned)) {
    stop("`adopt_unsigned` must be TRUE or FALSE.")
  }
  if (!file.exists(path)) stop("File does not exist: ", path)
  ensure_qs2()
  ensure_openssl()

  before <- file.info(path)[, c("size", "mtime")]
  parts <- split_omp_signature(read_file_bytes(path))
  result <- function(status, message) list(status = status, message = message)
  if (!is.null(parts$mac)) {
    if (mac_equal(omp_mac(parts$body, key), parts$mac)) {
      return(result("valid", "Already signed."))
    }
    if (is.null(old) || !mac_equal(omp_mac(parts$body, old), parts$mac)) {
      return(result("refused", "The signature matches no key this store knows."))
    }
    status <- "resigned"
  } else {
    if (!isTRUE(adopt_unsigned)) {
      return(result("refused", "The file carries no signature."))
    }
    problem <- tryCatch({
      project_from_envelope(qs2::qs_deserialize(parts$body), path, untrusted = TRUE)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem)) return(result("refused", problem))
    status <- "signed"
  }
  # A writer that replaced the file while it was being checked wins: its
  # version is newer, and what was checked is no longer what is there.
  now <- file.info(path)[, c("size", "mtime")]
  if (!identical(before, now)) {
    return(result("refused", "The file changed while it was being checked."))
  }
  write_signed_omp(parts$body, key, path)
  result(status, if (status == "signed") "Signed." else "Re-signed under the current key.")
}

# Accept anything with the same major version as the current schema.
schema_is_supported <- function(version) {
  current <- as.integer(strsplit(OMP_SCHEMA_VERSION, ".", fixed = TRUE)[[1L]][[1L]])
  parts <- strsplit(version, ".", fixed = TRUE)[[1L]]
  if (length(parts) == 0L) return(FALSE)
  major <- suppressWarnings(as.integer(parts[[1L]]))
  !is.na(major) && major <= current
}
