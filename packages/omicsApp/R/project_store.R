# Per-user project store.
#
# In the ShinyProxy deployment every user gets a private host directory
# bind-mounted at `/data` inside their container, and `OMICSAPP_DATA_DIR`
# points at it. Outside a container (developer laptop, `R CMD check`) we
# fall back to a directory under `tempdir()` so nothing here ever writes
# to a user's home by surprise.
#
# Two kinds of file live in the store:
#   * `<slug>.omp`   -- projects the user explicitly saved
#   * `_autosave.omp` -- the rolling snapshot written by app_server()
#
# The autosave slot is reserved: `project_slug()` refuses to produce it
# so an explicit save can never clobber the crash-recovery copy.

AUTOSAVE_SLUG <- "_autosave"

#' Resolve the per-user project directory
#'
#' Reads `OMICSAPP_DATA_DIR` and creates the directory if needed.
#'
#' @param create Whether to create the directory when it does not exist.
#'
#' @return Absolute path to the project directory.
#' @keywords internal
#' @noRd
omicsapp_data_dir <- function(create = TRUE) {
  dir <- Sys.getenv("OMICSAPP_DATA_DIR", "")
  if (!nzchar(dir)) {
    dir <- file.path(tempdir(), "omicsApp-projects")
  }
  if (isTRUE(create) && !dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir
}

#' Storage quota for the current user, in gigabytes
#'
#' Read from `OMICSAPP_QUOTA_GB`. A non-numeric or non-positive value
#' disables the check (returns `Inf`) rather than locking the user out
#' of their own data because of a typo in the deployment config.
#'
#' @return Numeric scalar; `Inf` when no quota is configured.
#' @keywords internal
#' @noRd
omicsapp_quota_gb <- function() {
  raw <- Sys.getenv("OMICSAPP_QUOTA_GB", "")
  if (!nzchar(raw)) return(Inf)
  value <- suppressWarnings(as.numeric(raw))
  if (is.na(value) || value <= 0) return(Inf)
  value
}

#' Bytes currently used by the project directory
#'
#' @param dir Project directory.
#'
#' @return Numeric scalar (bytes). `0` when the directory is absent.
#' @keywords internal
#' @noRd
data_dir_usage_bytes <- function(dir = omicsapp_data_dir()) {
  if (!dir.exists(dir)) return(0)
  files <- list.files(dir, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (length(files) == 0L) return(0)
  sizes <- file.info(files)$size
  sum(sizes, na.rm = TRUE)
}

#' Is the store at or over quota?
#'
#' @param dir Project directory.
#'
#' @return `TRUE` when a quota is configured and already reached.
#' @keywords internal
#' @noRd
quota_exceeded <- function(dir = omicsapp_data_dir()) {
  quota <- omicsapp_quota_gb()
  if (!is.finite(quota)) return(FALSE)
  (data_dir_usage_bytes(dir) / 1024^3) >= quota
}

#' Human-readable usage summary
#'
#' @param dir Project directory.
#'
#' @return A single string such as `"2.4 GB / 50 GB"`, or `"2.4 GB"`
#'   when no quota is configured.
#' @keywords internal
#' @noRd
# A size a person can read. The number matters here -- it is what tells
# someone whether deleting a project actually gave them room back.
format_bytes <- function(bytes) {
  if (!is.finite(bytes) || bytes <= 0) return("0 B")
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- min(length(units), 1L + floor(log(bytes, 1024)))
  sprintf(if (i == 1L) "%.0f %s" else "%.1f %s",
          bytes / 1024^(i - 1L), units[[i]])
}

usage_label <- function(dir = omicsapp_data_dir()) {
  used_gb <- data_dir_usage_bytes(dir) / 1024^3
  quota <- omicsapp_quota_gb()
  if (!is.finite(quota)) {
    return(sprintf("%.2f GB used", used_gb))
  }
  sprintf("%.2f GB / %g GB", used_gb, quota)
}

#' Is a string valid UTF-8?
#'
#' @param x A single string.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
is_valid_utf8 <- function(x) {
  y <- x
  Encoding(y) <- "UTF-8"
  !is.na(iconv(y, "UTF-8", "UTF-8"))
}

#' Truncate a string to a byte budget without splitting a character
#'
#' Filesystem name limits are counted in bytes (255 on ext4), but
#' `substr()` counts *characters* when the string is marked UTF-8 — so
#' an 80-character cap on CJK text yields a 240-byte filename. Cutting
#' the raw bytes instead can land mid-character, so we walk the tail
#' back until the result is valid UTF-8 again.
#'
#' @param x A single string.
#' @param max_bytes Byte budget.
#'
#' @return `x`, truncated to at most `max_bytes` bytes.
#' @keywords internal
#' @noRd
truncate_bytes <- function(x, max_bytes) {
  if (nchar(x, type = "bytes") <= max_bytes) return(x)
  enc <- Encoding(x)
  bytes <- utils::head(charToRaw(x), max_bytes)
  rebuild <- function(b) {
    out <- rawToChar(b)
    if (!identical(enc, "unknown")) Encoding(out) <- enc
    out
  }
  out <- rebuild(bytes)
  # At most 3 iterations: a UTF-8 sequence is 4 bytes at the widest.
  while (length(bytes) > 0L && !is_valid_utf8(out)) {
    bytes <- bytes[-length(bytes)]
    out <- rebuild(bytes)
  }
  out
}

#' Turn a user-supplied project name into a safe file stem
#'
#' Removes the characters that make a filename dangerous or unportable
#' rather than whitelisting the ones that are safe: a whitelist would
#' have to enumerate every script a user might name a project in, and a
#' CJK range in a PCRE class is locale-sensitive. The blacklist is pure
#' ASCII, so it behaves identically everywhere, and non-Latin names
#' survive untouched.
#'
#' Cleared: control characters, path separators, Windows-reserved
#' punctuation, and any `..` run. Leading dots and underscores are
#' stripped so the result can be neither a hidden file nor the reserved
#' autosave slot, and the stem is capped at 80 *bytes* — filesystem
#' name limits are counted in bytes, not characters.
#'
#' Every pattern below is spelled with explicit ASCII ranges rather than
#' POSIX classes like `[[:cntrl:]]`. UTF-8 is self-synchronising — a
#' multi-byte sequence never contains a byte below 0x80 — so an ASCII
#' class is safe even when R is matching byte-wise, which is what
#' happens under `LC_CTYPE=C`. A POSIX class is *not* safe there: it
#' shreds any non-Latin name into replacement characters, and container
#' base images do not reliably set a UTF-8 locale.
#'
#' @param name Raw name typed by the user.
#'
#' @return A slug string, or `NA_character_` when nothing usable remains.
#' @keywords internal
#' @noRd
project_slug <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    return(NA_character_)
  }
  slug <- gsub("^[\x01-\x20]+|[\x01-\x20]+$", "", name)
  slug <- gsub("[\x01-\x1f\x7f]", "", slug)
  slug <- gsub("[/\\\\]", "_", slug)
  slug <- gsub("[<>:\"|?*]", "_", slug)
  slug <- gsub("[.]{2,}", "_", slug)
  slug <- gsub("[ \t\r\n\f\v]+", "_", slug)
  slug <- gsub("_{2,}", "_", slug)
  slug <- gsub("^[._]+", "", slug)
  slug <- gsub("[._]+$", "", slug)
  slug <- truncate_bytes(slug, 80L)
  # Again after truncating: the cut can land on a dot, and `<stem>.` +
  # ".omp" gives `stem..omp` -- the double-dot shape stripped earlier
  # for safety, reintroduced by the length cap.
  slug <- gsub("[._]+$", "", slug)
  # Belt and braces against a future edit to the rules above
  # reintroducing traversal. Deliberately *not* `basename()`: that
  # converts to the native encoding and therefore throws outright on a
  # non-Latin name under LC_CTYPE=C, which is a locale container base
  # images really do ship with. Asserting fails closed and needs no
  # conversion.
  if (grepl("[/\\\\]", slug) || grepl("[.][.]", slug)) {
    return(NA_character_)
  }
  if (!nzchar(slug) || identical(slug, AUTOSAVE_SLUG)) {
    return(NA_character_)
  }
  slug
}

#' A path the operating system will take in any locale
#'
#' A name typed in the browser arrives marked UTF-8. R translates a
#' marked path to the native encoding before every file call, and in a
#' process whose locale is C that translation fails for anything
#' outside ASCII: "unable to translate ... to native encoding", from
#' save, from archive, from open. A container without a configured
#' locale is exactly that process. The filesystem, on the other hand,
#' takes bytes, and UTF-8 bytes are what every other tool on the host
#' writes -- so the mark is dropped and the bytes go through as they
#' are. `list_saved_projects()` puts the mark back on the way out.
#'
#' @param path A path, possibly marked UTF-8.
#' @return The same bytes, unmarked.
#' @keywords internal
#' @noRd
os_path <- function(path) {
  path <- enc2utf8(path)
  Encoding(path) <- "unknown"
  path
}

# The inverse: names read back from the filesystem are UTF-8 bytes,
# and saying so is what lets the JSON layer send them to the browser
# as characters rather than as <e8><9b><8b>.
as_utf8_name <- function(x) {
  Encoding(x) <- "UTF-8"
  x
}

#' Absolute path for a project slug inside the store
#'
#' @param slug Slug from `project_slug()`.
#' @param dir Project directory.
#'
#' @return Absolute file path ending in `.omp`.
#' @keywords internal
#' @noRd
project_path <- function(slug, dir = omicsapp_data_dir()) {
  os_path(file.path(dir, paste0(slug, ".omp")))
}

#' Path to the rolling autosave snapshot
#'
#' @param dir Project directory.
#'
#' @return Absolute file path.
#' @keywords internal
#' @noRd
autosave_path <- function(dir = omicsapp_data_dir()) {
  file.path(dir, paste0(AUTOSAVE_SLUG, ".omp"))
}

#' List the projects saved in the store
#'
#' The reserved autosave snapshot is excluded — it is surfaced through
#' its own "restore" affordance, not as a normal project.
#'
#' @param dir Project directory.
#'
#' @return A `data.frame` with columns `slug`, `path`, `size_mb`,
#'   `modified`, ordered most-recently-modified first. Zero rows when
#'   the store is empty.
#' @keywords internal
#' @noRd
list_saved_projects <- function(dir = omicsapp_data_dir()) {
  empty <- data.frame(slug = character(0), path = character(0),
                      size_mb = numeric(0),
                      modified = as.POSIXct(character(0)),
                      stringsAsFactors = FALSE)
  if (!dir.exists(dir)) return(empty)
  files <- list.files(dir, pattern = "\\.omp$", full.names = TRUE)
  files <- files[basename(files) != paste0(AUTOSAVE_SLUG, ".omp")]
  if (length(files) == 0L) return(empty)
  info <- file.info(files)
  # Marked before the regex runs: sub() on unmarked non-ASCII bytes in a
  # C locale warns that it cannot translate them, and it is right.
  names <- as_utf8_name(basename(files))
  out <- data.frame(
    slug     = sub("\\.omp$", "", names),
    path     = files,
    size_mb  = round(info$size / 1024^2, 2),
    modified = info$mtime,
    stringsAsFactors = FALSE
  )
  out[order(out$modified, decreasing = TRUE), , drop = FALSE]
}

#' Write a project into the store
#'
#' Wraps `omicsCore::save_project()` with the quota check and slug
#' handling. Quota is only enforced for *new* files: overwriting an
#' existing project cannot grow the store without bound, and refusing it
#' would leave a user who is over quota unable to tidy up.
#'
#' @param project An `omics_project`.
#' @param slug Slug from `project_slug()`.
#' @param dir Project directory.
#' @param overwrite Whether to replace an existing file.
#'
#' @return A list with `ok` (logical), `path`, and `message`.
#' @keywords internal
#' @noRd
store_save_project <- function(project, slug, dir = omicsapp_data_dir(),
                               overwrite = FALSE) {
  fail <- function(msg) list(ok = FALSE, path = NA_character_, message = msg)
  if (is.na(slug)) {
    return(fail("Please enter a project name using letters, digits, or spaces."))
  }
  if (!omicsCore::is_omics_project(project)) {
    return(fail("Nothing to save yet \u2014 import a file first."))
  }
  path <- project_path(slug, dir)
  exists_already <- file.exists(path)
  if (exists_already && !isTRUE(overwrite)) {
    return(fail(sprintf("'%s' already exists. Tick overwrite to replace it.", slug)))
  }
  if (!exists_already && quota_exceeded(dir)) {
    return(fail(sprintf(
      "Storage quota reached (%s). Delete a project before saving a new one.",
      usage_label(dir)
    )))
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  result <- tryCatch(
    {
      omicsCore::save_project(project, path, overwrite = TRUE)
      list(ok = TRUE, path = path,
           message = sprintf("Saved '%s'.", slug))
    },
    error = function(e) fail(paste0("Save failed: ", conditionMessage(e)))
  )
  result
}

#' Delete a saved project
#'
#' @param slug Slug to remove.
#' @param dir Project directory.
#'
#' @return A list with `ok` (logical) and `message`.
#' @keywords internal
#' @noRd
store_delete_project <- function(slug, dir = omicsapp_data_dir()) {
  if (is.na(slug) || !nzchar(slug)) {
    return(list(ok = FALSE, message = "No project selected."))
  }
  path <- project_path(slug, dir)
  if (!file.exists(path)) {
    return(list(ok = FALSE, message = sprintf("'%s' no longer exists.", slug)))
  }

  # Read before removing: the project is the only record of which
  # uploads belong to it.
  mine <- project_fingerprints(path)

  removed <- suppressWarnings(file.remove(path))
  if (!isTRUE(removed)) {
    return(list(ok = FALSE, message = sprintf("Could not delete '%s'.", slug)))
  }

  freed <- prune_orphan_uploads(mine, dir)
  msg <- if (freed$n > 0L) {
    sprintf("Deleted '%s' and %d archived upload(s), %s freed.",
            slug, freed$n, freed$label)
  } else {
    sprintf("Deleted '%s'.", slug)
  }
  list(ok = TRUE, message = msg)
}

# The fingerprints of every layer in a saved project, as the 12-character
# digests store_raw_upload() names its files with.
project_fingerprints <- function(path) {
  proj <- tryCatch(omicsCore::load_project(path), error = function(e) NULL)
  if (is.null(proj)) return(character(0))
  fps <- vapply(proj$experiments,
                function(e) e$source_fingerprint %||% NA_character_,
                character(1L))
  fps <- fps[!is.na(fps) & nzchar(fps)]
  unique(substr(sub(":.*$", "", fps), 1L, 12L))
}

# Archived uploads are what actually fill the quota -- a project is a few
# MB, the workbook behind it can be 170 -- so deleting a project without
# them frees nothing the user can see.
#
# Only the ones nothing else refers to. Two projects built from one file
# share its digest, and removing the archive because one of them went
# would quietly take it from the other.
prune_orphan_uploads <- function(digests, dir = omicsapp_data_dir()) {
  none <- list(n = 0L, label = "0 B")
  if (length(digests) == 0L) return(none)

  # Through list_saved_projects() rather than a glob: projects are `.omp`,
  # not `.rds`, and a pattern written from memory matched nothing --
  # which made every upload look orphaned and deleted a file another
  # project was still using.
  #
  # The autosave is deliberately included: it holds a real project, and
  # an upload it refers to is one a restore would need.
  others <- c(list_saved_projects(dir)$path,
              file.path(dir, paste0(AUTOSAVE_SLUG, ".omp")))
  others <- others[file.exists(others)]
  still_used <- unlist(lapply(others, project_fingerprints), use.names = FALSE)
  orphans <- setdiff(digests, still_used %||% character(0))
  if (length(orphans) == 0L) return(none)

  raw <- raw_dir(dir, create = FALSE)
  if (!dir.exists(raw)) return(none)
  files <- list.files(raw, full.names = TRUE)
  hit <- files[grepl(paste0("__(", paste(orphans, collapse = "|"), ")\\."),
                     basename(files))]
  if (length(hit) == 0L) return(none)

  bytes <- sum(file.size(hit), na.rm = TRUE)
  ok <- suppressWarnings(file.remove(hit))
  list(n = sum(ok), label = format_bytes(bytes))
}

#' Read a saved project back out of the store
#'
#' @param slug Slug to load.
#' @param dir Project directory.
#'
#' @return A list with `ok`, `project`, and `message`.
#' @keywords internal
#' @noRd
store_load_project <- function(slug, dir = omicsapp_data_dir()) {
  fail <- function(msg) list(ok = FALSE, project = NULL, message = msg)
  if (is.na(slug) || !nzchar(slug)) return(fail("No project selected."))
  path <- project_path(slug, dir)
  if (!file.exists(path)) {
    return(fail(sprintf("'%s' no longer exists.", slug)))
  }
  tryCatch(
    list(ok = TRUE, project = omicsCore::load_project(path),
         message = sprintf("Opened '%s'.", slug)),
    error = function(e) fail(paste0("Open failed: ", conditionMessage(e)))
  )
}

#' Write the rolling autosave snapshot
#'
#' Failures are swallowed deliberately: autosave is a convenience, and a
#' full disk or a read-only mount must never break the analysis the user
#' is in the middle of. The return value lets callers log if they care.
#'
#' @param project An `omics_project`.
#' @param dir Project directory.
#'
#' @return `TRUE` on success, `FALSE` otherwise.
#' @keywords internal
#' @noRd
store_autosave <- function(project, dir = omicsapp_data_dir()) {
  if (!omicsCore::is_omics_project(project)) return(FALSE)
  # The snapshot replaces itself in place, so it does not grow the
  # store; enforcing quota here would only strand a user's recovery
  # copy at exactly the moment they most need it.
  tryCatch(
    {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      }
      omicsCore::save_project(project, autosave_path(dir), overwrite = TRUE)
      TRUE
    },
    error = function(e) FALSE
  )
}

#' Read the autosave snapshot, if one exists
#'
#' @param dir Project directory.
#'
#' @return An `omics_project`, or `NULL` when absent or unreadable.
#' @keywords internal
#' @noRd
store_read_autosave <- function(dir = omicsapp_data_dir()) {
  path <- autosave_path(dir)
  if (!file.exists(path)) return(NULL)
  tryCatch(omicsCore::load_project(path), error = function(e) NULL)
}

#' Modification time of the autosave snapshot
#'
#' @param dir Project directory.
#'
#' @return A `POSIXct`, or `NULL` when there is no snapshot.
#' @keywords internal
#' @noRd
autosave_mtime <- function(dir = omicsapp_data_dir()) {
  path <- autosave_path(dir)
  if (!file.exists(path)) return(NULL)
  file.info(path)$mtime
}

#' Directory holding archived raw uploads
#'
#' @param dir Project directory.
#' @param create Whether to create the directory when it does not exist.
#'
#' @return Absolute path to the `raw/` sub-directory.
#' @keywords internal
#' @noRd
raw_dir <- function(dir = omicsapp_data_dir(), create = TRUE) {
  path <- file.path(dir, "raw")
  if (isTRUE(create) && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

#' Archive the file an import was parsed from
#'
#' `.omp` carries the parsed matrices, which is enough to keep working,
#' but not enough to answer "what did the file we started from look
#' like?" — a question that outlives the session, and that a methods
#' section eventually asks. Keeping the upload closes that gap.
#'
#' The stored name is derived from the fingerprint, so re-importing the
#' same file twice is idempotent rather than accumulating copies.
#'
#' Unlike the autosave snapshot this genuinely grows the store, so it is
#' subject to quota — and it fails soft: archiving is a convenience, and
#' losing it must never block an import the user is in the middle of.
#'
#' @param path Path to the uploaded file (Shiny's temp copy).
#' @param name Original file name, used for the readable part of the
#'   stored name.
#' @param fingerprint Fingerprint from `input_fingerprint()`.
#' @param dir Project directory.
#'
#' @return A list with `ok`, `path`, and `message`.
#' @keywords internal
#' @noRd
store_raw_upload <- function(path, name, fingerprint,
                             dir = omicsapp_data_dir()) {
  skip <- function(msg) list(ok = FALSE, path = NA_character_, message = msg)
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(skip("No uploaded file to archive."))
  }
  if (is.null(fingerprint) || is.na(fingerprint)) {
    return(skip("Upload has no fingerprint; not archived."))
  }
  digest <- substr(sub(":.*$", "", fingerprint), 1L, 12L)
  stem <- project_slug(tools::file_path_sans_ext(name %||% "upload"))
  if (is.na(stem)) stem <- "upload"
  ext <- tools::file_ext(name %||% "")
  target <- os_path(file.path(raw_dir(dir),
                              paste0(stem, "__", digest,
                                     if (nzchar(ext)) paste0(".", ext) else "")))
  if (file.exists(target)) {
    return(list(ok = TRUE, path = target,
                message = "Already archived."))
  }
  if (quota_exceeded(dir)) {
    return(skip(sprintf(
      "Storage quota reached (%s); the raw file was not archived.",
      usage_label(dir)
    )))
  }
  tryCatch(
    {
      ok <- file.copy(path, target, overwrite = FALSE)
      if (!isTRUE(ok)) return(skip("Could not archive the raw file."))
      list(ok = TRUE, path = target, message = "Raw file archived.")
    },
    error = function(e) skip(paste0("Could not archive the raw file: ",
                                    conditionMessage(e)))
  )
}

#' Snapshot a project to disk whenever it changes
#'
#' Every import and every completed analysis lands a new value on
#' `current_project`, so one observer covers both triggers.
#'
#' Writes on every change rather than debouncing. A debounce was tried
#' and removed: measured against the real reactive graph it saved
#' exactly one write, in the one scenario where a re-import replaces a
#' layer that already carries results, and the bytes that reached disk
#' were identical either way. What it cost was a window in which a
#' change existed only in memory — and this app is deployed behind
#' ShinyProxy, which stops containers with a signal that runs no R
#' handler. Trading ~0.1s of I/O for a window of silent data loss is
#' the wrong way round when the container is expected to die.
#'
#' Hangs off the project *state* rather than each analysis-completion
#' *event*. `current_project` is updated by the bundle-attach observer
#' in `app_server()`, so an observer firing on the same flush as a
#' finished bundle would read the project from before that bundle was
#' folded in, and persist a snapshot missing the very result that
#' triggered it. Reading the state that triggered you has no such
#' ordering hazard.
#'
#' Split out of `app_server()` so `writer` can be injected — the path
#' deciding whether a recycled container costs the user their work
#' should not be the untestable one.
#'
#' @param current_project A reactive (or `reactiveVal`) yielding the
#'   live `omics_project` or `NULL`.
#' @param writer Function called with the project. Injected rather than
#'   mocked so a test can count writes without depending on how the
#'   test runner resolves namespace bindings.
#'
#' @return Invisibly, the observer handle.
#' @keywords internal
#' @noRd
wire_autosave <- function(current_project, writer = store_autosave) {
  # Reads `current_project` without writing it, so this cannot
  # re-trigger itself.
  invisible(shiny::observe({
    proj <- current_project()
    if (is.null(proj)) return()
    # `store_autosave()` already swallows its own failures, but the
    # invariant belongs here too: an observer that throws is destroyed,
    # so one unlucky write would silently switch autosaving off for the
    # rest of the session -- the failure mode this whole path exists to
    # prevent, arrived at by a different road.
    tryCatch(writer(proj), error = function(e) NULL)
  }))
}
