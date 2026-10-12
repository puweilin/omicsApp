# The newest autosave snapshot in a store. Each session writes its own
# (`_autosave-<id>.omp`), so tests look for the newest rather than a
# fixed name; a path that does not exist when there is none.
autosave_file <- function(store) {
  files <- list.files(store, pattern = "^_autosave(-[A-Za-z0-9]+)?\\.omp$",
                      full.names = TRUE)
  if (!length(files)) return(file.path(store, "_autosave-none.omp"))
  files[order(file.info(files)$mtime, decreasing = TRUE)][[1L]]
}
