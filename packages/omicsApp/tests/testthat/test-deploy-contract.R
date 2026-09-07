# The deployment files have to agree with each other and with the
# packages, and nothing checked that they did. A third of this project's
# fix commits are deploy fixes, and most of them were the same failure:
# two files that had to agree quietly did not (the nginx port and the
# ShinyProxy port; the password floor in three places; the package
# lists in the Dockerfile and check_pins.R), or a setting was latent
# because no container had ever started. Each of those is a line here.
#
# These read the deploy/ directory at the repository root, so they run
# from the source tree and skip from an installed package.

deploy_root <- function() {
  dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    candidate <- file.path(dir, "deploy")
    if (file.exists(file.path(candidate, "docker", "Dockerfile"))) return(candidate)
    parent <- dirname(dir)
    if (identical(parent, dir)) return(NA_character_)
    dir <- parent
  }
}

skip_unless_deploy <- function() {
  root <- deploy_root()
  skip_if(is.na(root), "deploy/ directory not found; running from an installed package")
  root
}

read_deploy <- function(root, ...) readLines(file.path(root, ...), warn = FALSE)

# `key: value` from a YAML-ish file, ignoring comments. Good enough for
# the flat settings asserted here; a YAML parser would be another
# dependency for a test.
yaml_value <- function(lines, key) {
  hits <- grep(sprintf("^\\s*%s:\\s*", key), lines, value = TRUE)
  hits <- hits[!grepl("^\\s*#", hits)]
  if (length(hits) == 0L) return(NA_character_)
  trimws(sub(sprintf("^\\s*%s:\\s*", key), "", hits[[1L]]))
}

quoted_names <- function(lines) {
  unique(gsub("\"", "", unlist(regmatches(lines,
    gregexpr("\"[A-Za-z][A-Za-z0-9.]+\"", lines)))))
}

# The two package vectors in check_pins.R, read to their closing paren.
pinned_packages <- function(pins) {
  grab <- function(name) {
    i <- grep(sprintf("^%s <- c\\(", name), pins)
    if (length(i) != 1L) return(character(0))
    j <- i
    while (!grepl("\\)\\s*$", pins[j])) j <- j + 1L
    quoted_names(pins[i:j])
  }
  c(grab("CRAN_PKGS"), grab("BIOC_PKGS"))
}

# ---- the package lists --------------------------------------------------

test_that("check_pins.R and the Dockerfile install the same packages", {
  root <- skip_unless_deploy()
  pins <- read_deploy(root, "scripts", "check_pins.R")
  pinned <- pinned_packages(pins)
  expect_gt(length(pinned), 30L)
  docker <- quoted_names(read_deploy(root, "docker", "Dockerfile"))
  expect_length(setdiff(pinned, docker), 0L)
})

test_that("every package the image installs is one the packages declare", {
  root <- skip_unless_deploy()
  declared <- function(pkg) {
    d <- read.dcf(file.path(root, "..", "packages", pkg, "DESCRIPTION"))
    f <- function(field) {
      if (!field %in% colnames(d)) return(character(0))
      trimws(gsub("\\s*\\(.*?\\)", "", strsplit(d[1, field], ",")[[1]]))
    }
    c(f("Imports"), f("Suggests"))
  }
  wanted <- unique(c(declared("omicsCore"), declared("omicsApp")))
  pinned <- pinned_packages(read_deploy(root, "scripts", "check_pins.R"))
  # Everything pinned must be declared somewhere, or the image carries
  # a package nothing asks for.
  expect_length(setdiff(pinned, wanted), 0L)
  # The reverse is allowed only for the exclusions the Dockerfile
  # documents (ggpubr's dependency chain; tximport/GenomicFeatures are
  # left out of the image on purpose). Anything else missing would
  # silently remove a feature, since omicsCore gates on requireNamespace.
  documented_out <- c("ggpubr", "tximport", "GenomicFeatures")
  # Test-only packages the image has no use for. A package that belongs
  # in neither list is the finding this test exists for: it would be
  # missing from the image and, because omicsCore gates on
  # requireNamespace(), silently missing as a feature.
  dev_only <- c("testthat", "withr", "callr", "shinytest2", "chromote",
                "writexl", "later", "pkgload", "methods", "R", "omicsCore")
  expect_length(setdiff(wanted, c(pinned, documented_out, dev_only)), 0L)
})

# ---- the ports, which have to agree across three files -------------------

test_that("nginx forwards to the port ShinyProxy listens on", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  nginx <- read_deploy(root, "nginx", "omicsapp.conf.template")
  port <- yaml_value(sp, "port")
  expect_match(port, "^[0-9]+$")
  forwarded <- grep("proxy_pass http://127\\.0\\.0\\.1:[0-9]+;", nginx, value = TRUE)
  forwarded <- forwarded[!grepl("/auth/", forwarded)]
  expect_length(forwarded, 1L)
  expect_match(forwarded, paste0(":", port, ";"), fixed = TRUE)
})

test_that("ShinyProxy's two servers use different ports, both on localhost", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  main <- yaml_value(sp, "port")
  # management.server.port is the second `port:` key in the file
  ports <- grep("^\\s*port:\\s*[0-9]+", sp, value = TRUE)
  ports <- as.integer(sub(".*port:\\s*", "", ports))
  expect_gte(length(ports), 3L)   # proxy.port, spec container port, management
  management <- ports[[length(ports)]]
  expect_false(identical(as.integer(main), management))
  expect_identical(yaml_value(sp, "bind-address"), "127.0.0.1")
  addresses <- grep("^\\s*address:\\s*", sp, value = TRUE)
  expect_true(all(grepl("127.0.0.1", addresses, fixed = TRUE)))
})

test_that("the container port is the one launch() is told to use and the one exposed", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  docker <- read_deploy(root, "docker", "Dockerfile")
  cmd <- grep("^CMD", docker, value = TRUE)
  expect_length(cmd, 1L)
  expect_match(cmd, "host = '0.0.0.0'", fixed = TRUE)
  cmd_port <- sub(".*port = ([0-9]+).*", "\\1", cmd)
  exposed <- sub("^EXPOSE\\s+", "", grep("^EXPOSE", docker, value = TRUE))
  expect_identical(cmd_port, exposed)
  spec_port <- grep("^\\s+port:\\s*[0-9]+", sp, value = TRUE)
  expect_true(any(grepl(paste0(":\\s*", cmd_port, "$"), spec_port)))
})

# ---- the three faults that each cost a day --------------------------------

test_that("ShinyProxy is configured for a host install behind TLS", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  # Host-installed ShinyProxy cannot resolve Docker DNS names.
  expect_identical(yaml_value(sp, "internal-networking"), "false")
  # TLS is terminated by nginx; without these the OIDC redirect is http://.
  expect_identical(yaml_value(sp, "forward-headers-strategy"), "native")
  expect_identical(yaml_value(sp, "enforce-https-redirect-uri"), "true")
  expect_identical(yaml_value(sp, "secure-cookies"), "true")
  # The one claim Keycloak guarantees never changes names the volume.
  expect_identical(yaml_value(sp, "username-attribute"), "sub")
  expect_identical(yaml_value(sp, "roles-claim"), "groups")
  expect_identical(yaml_value(sp, "authentication"), "openid")
})

test_that("the servlet session outlives the heartbeat timeout", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  heartbeat_ms <- as.numeric(yaml_value(sp, "heartbeat-timeout"))
  timeout <- yaml_value(sp, "timeout")
  unit <- sub("^[0-9]+", "", timeout)
  n <- as.numeric(sub("[a-z]+$", "", timeout))
  session_ms <- n * switch(unit, h = 3600e3, m = 60e3, s = 1e3, 1)
  expect_gte(session_ms, heartbeat_ms)
})

test_that("nginx carries the WebSocket and the forwarded scheme, and no http2", {
  root <- skip_unless_deploy()
  nginx <- read_deploy(root, "nginx", "omicsapp.conf.template")
  code <- nginx[!grepl("^\\s*#", nginx)]
  expect_false(any(grepl("http2", code)))
  expect_true(any(grepl("proxy_set_header Upgrade\\s+\\$http_upgrade", code)))
  expect_true(any(grepl("proxy_set_header Connection\\s+\"upgrade\"", code)))
  expect_gte(sum(grepl("X-Forwarded-Proto\\s+\\$scheme", code)), 2L)
  expect_true(any(grepl("proxy_read_timeout", code)))
})

test_that("nginx accepts uploads at least as large as the app does", {
  root <- skip_unless_deploy()
  nginx <- read_deploy(root, "nginx", "omicsapp.conf.template")
  line <- grep("client_max_body_size", nginx, value = TRUE)
  expect_length(line, 1L)
  nginx_mb <- as.numeric(sub(".*client_max_body_size\\s+([0-9]+)M;.*", "\\1", line))
  app_mb <- formals(launch)$max_upload_mb
  expect_gte(nginx_mb, app_mb)
})

# ---- the password floor, which drifted across three files -----------------

test_that("the password policy is the same number everywhere it is written", {
  root <- skip_unless_deploy()
  realm <- read_deploy(root, "keycloak", "omicsapp-realm.json.template")
  policy <- grep("\"passwordPolicy\"", realm, value = TRUE)
  expect_length(policy, 1L)
  floor <- as.integer(sub(".*length\\(([0-9]+)\\).*", "\\1", policy))

  script <- read_deploy(root, "scripts", "add_user.sh")
  checks <- grep("len\\(SHARED\\) < [0-9]+", script, value = TRUE)
  expect_gte(length(checks), 1L)
  expect_true(all(as.integer(sub(".*< ([0-9]+).*", "\\1", checks)) == floor))

  readme <- read_deploy(root, "keycloak", "README.md")
  documented <- grep("`passwordPolicy`", readme, value = TRUE)
  expect_true(any(grepl(sprintf("length\\(%d\\)", floor), documented)))
})

# ---- the gene-set cache the image bakes in --------------------------------

test_that("the prewarm script builds every database omicsCore knows, for both organisms", {
  root <- skip_unless_deploy()
  prewarm <- read_deploy(root, "docker", "prewarm_genesets.R")
  block <- prewarm[grep("^COLLECTIONS <- list\\(", prewarm):length(prewarm)]
  block <- block[seq_len(grep("^\\)", block)[[1L]])]
  baked <- sub("^\\s*([a-z_]+)\\s*=.*", "\\1", grep("^\\s*[a-z_]+\\s*= list\\(", block, value = TRUE))
  expect_setequal(baked, names(omicsCore:::DB_MSIGDBR_MAP))

  organisms <- grep("^ORGANISMS <- c\\(", prewarm, value = TRUE)
  n_org <- lengths(regmatches(organisms, gregexpr("\"[^\"]+\"", organisms)))
  # The README tells the operator what count to expect after the build.
  readme <- read_deploy(root, "README.md")
  expected <- as.integer(sub(".*# expect ([0-9]+).*", "\\1",
                             grep("^# expect [0-9]+", readme, value = TRUE)[[1L]]))
  expect_identical(length(baked) * n_org, expected)

  docker <- read_deploy(root, "docker", "Dockerfile")
  cache_env <- grep("^ENV OMICSCORE_GENESET_CACHE=", docker, value = TRUE)
  expect_length(cache_env, 1L)
  default_dir <- sub('.*"OMICSCORE_GENESET_CACHE", "([^"]+)".*', "\\1",
                     grep('Sys.getenv("OMICSCORE_GENESET_CACHE"', prewarm,
                          value = TRUE, fixed = TRUE))
  expect_identical(sub("^ENV OMICSCORE_GENESET_CACHE=", "", cache_env), default_dir)
})

# ---- the storage layout, which three files name ----------------------------

test_that("the user store is the same directory in ShinyProxy, add_user.sh and the backup", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  volume <- grep("#\\{proxy.userId\\}:/data", sp, value = TRUE)
  volume <- volume[!grepl("^\\s*#", volume)]
  expect_length(volume, 1L)
  users_root <- sub('.*"([^"]+)/#\\{proxy.userId\\}:/data".*', "\\1", volume)

  script <- read_deploy(root, "scripts", "add_user.sh")
  default_root <- grep("^USERS_ROOT=", script, value = TRUE)
  expect_match(default_root, users_root, fixed = TRUE)

  cron <- read_deploy(root, "cron", "omicsapp-backup")
  expect_true(any(grepl(paste0(users_root, "/"), cron, fixed = TRUE)))

  # Inside the container the store is /data, and the app reads it from
  # this variable -- so the image must set it, or every project lands
  # in tempdir() and vanishes with the container.
  expect_true(any(grepl("OMICSAPP_QUOTA_GB", sp)))
  docker <- read_deploy(root, "docker", "Dockerfile")
  expect_true(any(grepl("OMICSAPP_DATA_DIR=/data", docker)))
})

test_that("the image runs in a UTF-8 locale", {
  root <- skip_unless_deploy()
  docker <- read_deploy(root, "docker", "Dockerfile")
  code <- docker[!grepl("^\\s*#", docker)]
  expect_true(any(grepl("LANG=[A-Za-z_]+\\.UTF-8", code)))
  expect_true(any(grepl("LC_ALL=[A-Za-z_]+\\.UTF-8", code)))
})

test_that("the CRAN snapshot is a date the check_pins script can read", {
  root <- skip_unless_deploy()
  docker <- read_deploy(root, "docker", "Dockerfile")
  snap <- sub("^ARG CRAN_SNAPSHOT=", "", grep("^ARG CRAN_SNAPSHOT=", docker, value = TRUE))
  expect_length(snap, 1L)
  expect_false(is.na(as.Date(snap, format = "%Y-%m-%d")))
})

# ---- the server's address, which four files carry ---------------------------
#
# It used to be written into all four by hand, and the day it changes
# they change together or not at all. Now it lives in host.env and
# render.sh writes the four. These check that no tracked deploy file
# names a site, that every carrier is a template with the token, and
# that rendering fills every token and leaves only the client secret.

rendered_outputs <- c("nginx/omicsapp.conf", "keycloak/docker-compose.yml",
                      "keycloak/omicsapp-realm.json", "shinyproxy/application.yml")

render_script <- function(root) file.path(root, "scripts", "render.sh")

skip_unless_bash <- function() {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("bash")), "bash not available")
}

test_that("no tracked deploy file names a site-specific address", {
  root <- skip_unless_deploy()
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  rel <- sub(paste0("^", root, "/"), "", files)
  # The rendered files and host.env are gitignored and may exist locally.
  files <- files[!rel %in% c(rendered_outputs, "host.env")]
  hits <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    i <- grep("192\\.168\\.[0-9]+\\.[0-9]+|(^|[^0-9.])10\\.[0-9]+\\.[0-9]+\\.[0-9]+", lines)
    if (length(i)) paste0(sub(paste0("^", root, "/"), "", f), ":", i) else character()
  }))
  expect_length(hits, 0L)
})

test_that("each carrier of the address is a template, rendered, ignored and kept from rsync", {
  root <- skip_unless_deploy()
  render <- read_deploy(root, "scripts", "render.sh")
  excl <- read_deploy(root, "rsync.exclude")
  ignore <- readLines(file.path(root, "..", ".gitignore"), warn = FALSE)
  for (out in rendered_outputs) {
    tmpl <- read_deploy(root, paste0(out, ".template"))
    expect_true(any(grepl("@OMICSAPP_HOST@", tmpl, fixed = TRUE)), info = out)
    expect_true(any(grepl(out, render, fixed = TRUE)), info = out)
    expect_true(paste0("deploy/", out) %in% excl, info = out)
    expect_true(paste0("deploy/", out) %in% ignore, info = out)
  }
  expect_true("deploy/host.env" %in% excl)
  expect_true("deploy/host.env" %in% ignore)
  expect_true("deploy/keycloak/.env" %in% excl)
  expect_true(any(grepl("^OMICSAPP_HOST=REPLACE_ME$", read_deploy(root, "host.env.template"))))
})

test_that("render.sh fills every token and leaves only the client secret to fill", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  out <- withr::local_tempdir()
  status <- system2("bash", c(render_script(root), "--host", "203.0.113.7", "--out", out),
                    stdout = FALSE, stderr = FALSE)
  expect_identical(status, 0L)
  for (f in rendered_outputs) {
    lines <- readLines(file.path(out, f), warn = FALSE)
    expect_false(any(grepl("@OMICSAPP_", lines, fixed = TRUE)), info = f)
    expect_true(any(grepl("203.0.113.7", lines, fixed = TRUE)), info = f)
  }
  nginx <- readLines(file.path(out, "nginx", "omicsapp.conf"), warn = FALSE)
  expect_length(grep("^\\s*server_name .*203\\.0\\.113\\.7;", nginx), 2L)
  expect_true(any(grepl("subjectAltName=IP:203.0.113.7,", nginx, fixed = TRUE)))
  compose <- readLines(file.path(out, "keycloak", "docker-compose.yml"), warn = FALSE)
  expect_identical(yaml_value(compose, "KC_HOSTNAME"), "https://203.0.113.7/auth")
  sp <- readLines(file.path(out, "shinyproxy", "application.yml"), warn = FALSE)
  expect_identical(yaml_value(sp, "auth-url"),
                   "https://203.0.113.7/auth/realms/omicsapp/protocol/openid-connect/auth")
  expect_identical(yaml_value(sp, "client-secret"), "REPLACE_ME")

  skip_if_not_installed("jsonlite")
  realm <- jsonlite::fromJSON(file.path(out, "keycloak", "omicsapp-realm.json"),
                              simplifyVector = FALSE)
  client <- realm$clients[[1L]]
  expect_identical(client$redirectUris[[1L]], "https://203.0.113.7/login/oauth2/code/shinyproxy")
  expect_identical(client$webOrigins[[1L]], "https://203.0.113.7")
  expect_identical(client$attributes$post.logout.redirect.uris, "https://203.0.113.7/*")
})

test_that("a DNS name gets a DNS subjectAltName, because Java matches the entry type", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  out <- withr::local_tempdir()
  status <- system2("bash", c(render_script(root), "--host", "omics.example.org", "--out", out),
                    stdout = FALSE, stderr = FALSE)
  expect_identical(status, 0L)
  nginx <- readLines(file.path(out, "nginx", "omicsapp.conf"), warn = FALSE)
  expect_true(any(grepl("subjectAltName=DNS:omics.example.org,", nginx, fixed = TRUE)))
  expect_length(grep("^\\s*server_name .*omics\\.example\\.org;", nginx), 2L)
})

test_that("render.sh refuses an address with a scheme, port or path, and the placeholder", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  out <- withr::local_tempdir()
  for (bad in c("https://10.0.0.1", "host:8443", "host/path", "REPLACE_ME")) {
    status <- system2("bash", c(render_script(root), "--host", bad, "--out", out),
                      stdout = FALSE, stderr = FALSE)
    expect_false(identical(status, 0L), info = bad)
  }
  expect_length(list.files(out, recursive = TRUE), 0L)
})

# --- the backup that could destroy the previous one ----------------------

# The crontab used to hold the whole thing:
#
#   docker exec keycloak-db pg_dump ... | gzip > /backup/.../keycloak-$(date +%u).sql.gz
#
# Two failures in one line. Cron's /bin/sh has no pipefail, so the status
# is gzip's, and gzip succeeds at compressing nothing. And `>` truncates
# the target before the pipeline runs, so a night when pg_dump could not
# start destroyed the backup from a week ago and left a 20-byte file that
# `gzip -t` calls valid.

backup_script <- function() {
  path <- file.path(deploy_root(), "scripts", "backup-keycloak.sh")
  skip_if_not(file.exists(path), "backup script not present")
  path
}

# A stand-in `docker` on PATH, so the script can be exercised without one.
fake_docker <- function(behaviour, env = parent.frame()) {
  bin <- withr::local_tempdir(.local_envir = env)
  writeLines(c("#!/bin/bash", behaviour), file.path(bin, "docker"))
  Sys.chmod(file.path(bin, "docker"), "0755")
  bin
}

run_backup <- function(bin, dest_dir) {
  withr::with_envvar(
    c(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
      BACKUP_DIR = dest_dir),
    suppressWarnings(system2("bash", shQuote(backup_script()),
                             stdout = TRUE, stderr = TRUE)))
}

test_that("the crontab no longer pipes into the rotation target", {
  cron <- readLines(file.path(deploy_root(), "cron", "omicsapp-backup"))
  active <- grep("^[^#]", cron, value = TRUE)
  expect_false(any(grepl("pg_dump.*\\|.*gzip.*>", active)),
               info = "the pipeline is back in the crontab, where /bin/sh has no pipefail")
  expect_true(any(grepl("backup-keycloak.sh", active)))
})

test_that("a good dump is written and verified", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  bin <- fake_docker('printf "PostgreSQL database dump\\n%.0s" {1..500}')
  dest <- withr::local_tempdir()
  out <- run_backup(bin, dest)
  expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(out, collapse = "\n"))
  written <- list.files(dest, pattern = "^keycloak-.*\\.sql\\.gz$", full.names = TRUE)
  expect_length(written, 1L)
  expect_gt(file.info(written)$size, 0)
  # And nothing is left beside it.
  expect_length(list.files(dest), 1L)
})

test_that("a failed dump leaves the previous backup byte-for-byte intact", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  dest <- withr::local_tempdir()

  ok_bin <- fake_docker('printf "PostgreSQL database dump\\n%.0s" {1..500}')
  expect_equal(attr(run_backup(ok_bin, dest), "status") %||% 0L, 0L)
  written <- list.files(dest, full.names = TRUE)[[1L]]
  before <- tools::md5sum(written)

  # pg_dump cannot start. The old line would have truncated the target
  # first and then reported success.
  bad_bin <- fake_docker('echo "could not connect to server" >&2; exit 1')
  out <- run_backup(bad_bin, dest)
  expect_gt(attr(out, "status") %||% 0L, 0L)
  expect_match(paste(out, collapse = " "), "previous")

  expect_true(file.exists(written))
  expect_identical(tools::md5sum(written), before)
  expect_length(list.files(dest, all.files = TRUE, no.. = TRUE), 1L)
})

test_that("an empty but valid gzip is refused", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  # Exactly what the old pipeline produced on a bad night: 20 bytes that
  # `gzip -t` accepts. Size alone is a coarse screen, so the script also
  # checks that what it wrote looks like pg_dump output.
  dest <- withr::local_tempdir()
  bin <- fake_docker('exit 0')       # succeeds, prints nothing
  out <- run_backup(bin, dest)
  expect_gt(attr(out, "status") %||% 0L, 0L)
  expect_length(list.files(dest, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("output that is not a database dump is refused", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  dest <- withr::local_tempdir()
  # Large enough to pass a size floor, and not a backup.
  bin <- fake_docker('printf "FATAL: password authentication failed\\n%.0s" {1..500}')
  out <- run_backup(bin, dest)
  expect_gt(attr(out, "status") %||% 0L, 0L)
  expect_match(paste(out, collapse = " "), "does not look like pg_dump")
  expect_length(list.files(dest, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("a large valid dump succeeds, which the small ones could not show", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  # The check used to stop reading at the header, so gunzip took SIGPIPE on
  # anything past a pipe buffer and `pipefail` failed every real backup
  # while passing every tiny one -- exactly backwards, and invisible to a
  # test whose fixture fits in the buffer.
  dest <- withr::local_tempdir()
  bin <- fake_docker('printf "PostgreSQL database dump\\n"; head -c 12000000 /dev/zero | tr "\\0" "x"')
  out <- run_backup(bin, dest)
  expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(out, collapse = "\n"))
  written <- list.files(dest, full.names = TRUE)
  expect_length(written, 1L)
  expect_gt(file.info(written)$size, 1000)
})

test_that("a large dump that is not a dump is still refused", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  # Large enough to clear the floor and past the buffer, so this also
  # shows the header check still runs to completion.
  dest <- withr::local_tempdir()
  bin <- fake_docker('head -c 12000000 /dev/zero | tr "\\0" "x"')
  out <- run_backup(bin, dest)
  expect_gt(attr(out, "status") %||% 0L, 0L)
  expect_match(paste(out, collapse = " "), "does not look like pg_dump")
  expect_length(list.files(dest, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("a compression failure does not overwrite the previous backup", {
  skip_on_os("windows")
  skip_if_not_installed("withr")
  dest <- withr::local_tempdir()
  ok_bin <- fake_docker('printf "PostgreSQL database dump\\n%.0s" {1..500}')
  expect_equal(attr(run_backup(ok_bin, dest), "status") %||% 0L, 0L)
  kept <- list.files(dest, full.names = TRUE)[[1L]]
  before <- tools::md5sum(kept)

  # gzip itself fails: a full disk on the backup volume looks like this.
  bin <- withr::local_tempdir()
  writeLines(c("#!/bin/bash",
               'printf "PostgreSQL database dump\\n%.0s" {1..500}'),
             file.path(bin, "docker"))
  writeLines(c("#!/bin/bash", "exit 1"), file.path(bin, "gzip"))
  Sys.chmod(file.path(bin, "docker"), "0755")
  Sys.chmod(file.path(bin, "gzip"), "0755")
  out <- run_backup(bin, dest)
  expect_gt(attr(out, "status") %||% 0L, 0L)
  expect_identical(tools::md5sum(kept), before)
})
