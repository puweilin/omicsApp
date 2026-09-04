# Running omicsCore in another R process.
#
# Two suites need it: the concurrent-writer test, which stages two
# processes on one file, and the locale test, which runs the package
# under LC_ALL=C the way a container without a configured locale does.
# Both have to load the *same* omicsCore this process is running --
# the source tree under load_all, or the library the package was loaded
# from when installed. Named explicitly, because R CMD check installs
# into a temporary library and a child that just says library(omicsCore)
# may pick up an older copy from the user library instead.

omicscore_origin <- function() {
  root <- system.file(package = "omicsCore")
  if (identical(basename(root), "inst")) root <- dirname(root)
  # An installed package has an R/ directory too -- it holds the
  # lazy-load database -- so the source tree is recognised by its .R
  # files, not by the directory.
  has_sources <- dir.exists(file.path(root, "R")) &&
    length(list.files(file.path(root, "R"), pattern = "\\.[Rr]$")) > 0L
  if (file.exists(file.path(root, "DESCRIPTION")) && has_sources &&
      requireNamespace("pkgload", quietly = TRUE)) {
    return(list(kind = "source", path = normalizePath(root)))
  }
  list(kind = "installed", path = dirname(normalizePath(root)))
}

# Code every child runs first. Returned as an expression-free function
# so callr can serialise it without dragging a test environment along.
load_omicscore_in_child <- function(origin) {
  if (identical(origin$kind, "source")) {
    pkgload::load_all(origin$path, quiet = TRUE)
  } else {
    library(omicsCore, lib.loc = c(origin$path, .libPaths()))
  }
  invisible(TRUE)
}

# Run `fn(origin, ...)` in a fresh R process. `env` is added to callr's
# safe defaults; pass LC_ALL and friends to change the locale.
in_child_process <- function(fn, ..., env = character(0), timeout = 300) {
  callr::r(fn, args = list(omicscore_origin(), ...),
           env = c(callr::rcmd_safe_env(), env), timeout = timeout)
}
