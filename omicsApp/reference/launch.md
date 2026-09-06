# Launch the omicsApp Shiny application

Starts the multi-omics analysis web interface.

## Usage

``` r
launch(
  port = NULL,
  host = "127.0.0.1",
  project = NULL,
  launch.browser = interactive(),
  workers = 2L,
  max_upload_mb = 500,
  ...
)
```

## Arguments

- port:

  Integer port number. If `NULL` (default), Shiny chooses a free port.

- host:

  Host address to bind to. Default `"127.0.0.1"` (local only). Use
  `"0.0.0.0"` when running inside a container behind a proxy.

- project:

  Optional path to an existing `.omp` project file to open on startup.

- launch.browser:

  Logical. If `TRUE` (default), open the app in the user's default
  browser. Set to `FALSE` for headless server environments.

- workers:

  Integer number of background workers for asynchronous analyses. `0`
  keeps the current `future` plan untouched (useful in tests, where
  `sequential` is what we want).

- max_upload_mb:

  Maximum accepted upload size in megabytes. Shiny's built-in default is
  5 MB, which is below the size of a typical omics expression workbook.

- ...:

  Additional arguments passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Invisible `NULL`. Called for its side effect (starts the Shiny server).

## Running as a multi-user server

Two settings matter when omicsApp is served rather than run on a laptop,
and both default to laptop-friendly values:

- `host = "0.0.0.0"` binds all interfaces so a reverse proxy (or a
  container port mapping) can reach the app. Keep the default
  `"127.0.0.1"` for local use.

- `workers` controls how many background R processes long analyses are
  offloaded to. Without this,
  [`future::future()`](https://future.futureverse.org/reference/future.html)
  in `run_async()` falls back to the `sequential` plan and every
  computation blocks the Shiny event loop — one user running DESeq2
  freezes every other session sharing the process.

On Linux (including inside a container) `workers > 0` selects
[`future::multicore`](https://future.futureverse.org/reference/multicore.html),
which forks and therefore shares the parent's expression matrices
copy-on-write. Where forking is unavailable or unsafe — Windows, RStudio
— it falls back to
[`future::multisession`](https://future.futureverse.org/reference/multisession.html),
which serialises captured data to each worker instead.
