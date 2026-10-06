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

  backup <- read_deploy(root, "scripts", "backup.sh")
  expect_true(any(grepl(sprintf('USERS_DIR="${USERS_DIR:-%s}"', users_root), backup, fixed = TRUE)))
  env_tpl <- read_deploy(root, "backup.env.template")
  expect_true(any(grepl(paste0("USERS_DIR=", users_root), env_tpl, fixed = TRUE)))

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

# ---- backups: versioned, verified, off-host, alerting ----------------------

test_that("the cron file runs the backup nightly and the restore drill weekly", {
  root <- skip_unless_deploy()
  cron <- read_deploy(root, "cron", "omicsapp-backup")
  jobs <- cron[!grepl("^\\s*(#|$)", cron) & grepl("root", cron)]
  expect_true(any(grepl("deploy/scripts/backup.sh", jobs, fixed = TRUE)))
  expect_true(any(grepl("deploy/scripts/restore_check.sh", jobs, fixed = TRUE)))
  # The pipe that replaced a good dump with an empty one is gone.
  expect_false(any(grepl("pg_dump", jobs, fixed = TRUE)))
})

test_that("the backup script fails loudly and keeps history", {
  root <- skip_unless_deploy()
  sh <- read_deploy(root, "scripts", "backup.sh")
  code <- sh[!grepl("^\\s*#", sh)]
  expect_true(any(grepl("set -Eeuo pipefail", code, fixed = TRUE)))
  expect_true(any(grepl("--link-dest", code, fixed = TRUE)))
  expect_true(any(grepl("BACKUP_REMOTE", code, fixed = TRUE)))
  expect_true(any(grepl("trap 'on_error", code, fixed = TRUE)))
  expect_true(any(grepl("MANIFEST.sha256", code, fixed = TRUE)))
})

test_that("backup.sh and restore_check.sh work end to end on a scratch layout", {
  root <- skip_unless_deploy()
  skip_on_os("windows")
  for (tool in c("bash", "rsync", "flock", "sha256sum", "gzip")) {
    skip_if(!nzchar(Sys.which(tool)), paste(tool, "not available"))
  }
  t <- withr::local_tempdir()
  dir.create(file.path(t, "users", "u1"), recursive = TRUE)
  dir.create(file.path(t, "kc"))
  dir.create(file.path(t, "remote"))
  writeLines("project", file.path(t, "users", "u1", "a.omp"))
  writeLines("16", file.path(t, "kc", "PG_VERSION"))
  dir.create(file.path(t, "logs", "shinyproxy", "containers"), recursive = TRUE)
  writeLines("R said something", file.path(t, "logs", "shinyproxy", "containers", "s1.log"))
  writeLines("backup ran", file.path(t, "logs", "omicsapp-backup.log"))
  docker <- file.path(t, "docker")
  writeLines(c("#!/bin/sh",
               "case \"$1\" in",
               "  exec) [ -n \"$FAKE_DUMP_FAIL\" ] && exit 1; head -c 3000 /dev/urandom | base64 ;;",
               "  image) echo 'sha256:test omicsapp:1.0' ;;",
               "esac"), docker)
  Sys.chmod(docker, "755")
  env <- c(OMICSAPP_BACKUP_CONF = "/nonexistent", USERS_DIR = file.path(t, "users"),
           KEYCLOAK_DB_DIR = file.path(t, "kc"), GENESETS_DIR = file.path(t, "none"),
           REPO_DIR = t, BACKUP_ROOT = file.path(t, "backup"), DOCKER = docker,
           ALERT_LOG = file.path(t, "alerts"), EXTRA_PATHS = "",
           LOG_PATHS = paste(file.path(t, "logs", "shinyproxy"),
                             file.path(t, "logs", "omicsapp-backup.log"),
                             file.path(t, "logs", "absent")),
           SP_CONF = "/nonexistent",
           BACKUP_REMOTE = file.path(t, "remote"), SKIP_DOCKER_CHECKS = "1")
  run <- function(script, extra = character(0)) {
    withr::with_envvar(c(env, extra),
      system2("bash", file.path(root, "scripts", script), stdout = FALSE, stderr = FALSE))
  }
  expect_identical(run("backup.sh"), 0L)
  Sys.sleep(1.1)
  writeLines("more", file.path(t, "users", "u1", "b.omp"))
  expect_identical(run("backup.sh"), 0L)
  snaps <- list.files(file.path(t, "backup", "snapshots"))
  expect_length(snaps, 2L)
  # Unchanged files are shared between snapshots, not copied.
  links <- system2("stat", c("-c", "%h", shQuote(file.path(t, "backup", "snapshots",
                                                      snaps[[1L]], "users", "u1", "a.omp"))),
                   stdout = TRUE)
  expect_identical(as.integer(links), 2L)
  expect_setequal(list.files(file.path(t, "remote", "snapshots")), snaps)
  # The logs travel with the snapshot, under their own absolute paths.
  logs <- file.path(t, "backup", "latest", "logs", t, "logs")
  expect_true(file.exists(file.path(logs, "shinyproxy", "containers", "s1.log")))
  expect_true(file.exists(file.path(logs, "omicsapp-backup.log")))
  expect_identical(run("restore_check.sh"), 0L)

  # A failed dump alerts, exits non-zero, and leaves the snapshots alone.
  expect_false(run("backup.sh", c(FAKE_DUMP_FAIL = "1")) == 0L)
  expect_length(list.files(file.path(t, "backup", "snapshots")), 2L)
  expect_true(any(grepl("failed", readLines(file.path(t, "alerts")))))

  # A damaged file is caught by the drill.
  cat("x", file = file.path(t, "backup", "latest", "users", "u1", "b.omp"), append = TRUE)
  expect_false(run("restore_check.sh") == 0L)
})

# ---- CI builds and tests the image the server runs ----------------------------

test_that("CI builds the production Dockerfile and tests inside it", {
  root <- skip_unless_deploy()
  wf <- file.path(dirname(root), ".github", "workflows", "production-image.yaml")
  skip_if(!file.exists(wf), "workflows not present")
  lines <- readLines(wf, warn = FALSE)
  expect_true(any(grepl("file: deploy/docker/Dockerfile", lines, fixed = TRUE)))
  expect_true(any(grepl("pull_request", lines, fixed = TRUE)))
  expect_true(any(grepl('test_local(\\"packages/omicsCore\\"', lines, fixed = TRUE)))
  expect_true(any(grepl('test_local(\\"packages/omicsApp\\"', lines, fixed = TRUE)))
})

# ---- admin endpoints: from the allow-list only --------------------------------

uncommented <- function(lines) lines[!grepl("^\\s*#", lines)]

test_that("nginx refuses the admin endpoints to anyone not on the allow-list", {
  root <- skip_unless_deploy()
  nginx <- uncommented(read_deploy(root, "nginx", "omicsapp.conf.template"))
  # The allow-list is rendered into a geo block; the default is "not allowed".
  geo <- grep("^geo \\$omicsapp_admin_client", nginx)
  expect_length(geo, 1L)
  expect_match(nginx[geo + 1L], "default 0;", fixed = TRUE)
  expect_match(nginx[geo + 2L], "@OMICSAPP_ADMIN_ALLOW@", fixed = TRUE)
  # Every administrative path is in the map ...
  map <- nginx[grep("^map \\$uri \\$omicsapp_admin_path", nginx):length(nginx)]
  map <- map[seq_len(grep("^\\}", map)[[1L]])]
  for (path in c("/auth/admin", "/auth/realms/master", "/auth/(health|metrics)",
                 "/admin", "/actuator")) {
    expect_true(any(grepl(paste0("~*^", path, "(/.*)?$"), map, fixed = TRUE)), info = path)
  }
  # ... and the users' own account page is not.
  expect_false(any(grepl("realms/omicsapp", map, fixed = TRUE)))
  expect_false(any(grepl("realms/(", map, fixed = TRUE)))
  # Matched on the normalised $uri, which is what proxy_pass forwards.
  expect_false(any(grepl("\\$request_uri", map)))
  # Refused at server level, before either location is chosen, in the
  # HTTPS server (the HTTP one only redirects).
  https <- nginx[grep("listen 443 ssl;", nginx, fixed = TRUE):length(nginx)]
  guard <- grep("if \\(\\$omicsapp_admin_denied\\) \\{", https)
  expect_length(guard, 1L)
  expect_match(https[guard + 1L], "return 403;", fixed = TRUE)
  expect_lt(guard, grep("location /auth/", https, fixed = TRUE)[[1L]])
})

test_that("render.sh writes the allow-list: the server itself by default, or what is set", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  geo_line <- function(out) {
    nginx <- readLines(file.path(out, "nginx", "omicsapp.conf"), warn = FALSE)
    trimws(nginx[grep("^geo \\$omicsapp_admin_client", nginx) + 2L])
  }
  render <- function(out, ...) {
    system2("bash", c(render_script(root), "--host", "203.0.113.7", "--out", out, ...),
            stdout = FALSE, stderr = FALSE)
  }
  withr::local_envvar(ADMIN_ALLOW_CIDR = NA)
  out <- withr::local_tempdir()
  expect_identical(render(out), 0L)
  expect_identical(geo_line(out), "127.0.0.1 1; ::1 1; 203.0.113.7 1;")

  out2 <- withr::local_tempdir()
  expect_identical(render(out2, "--admin-allow", shQuote("198.51.100.0/24, 2001:db8::/32")), 0L)
  expect_identical(geo_line(out2), "198.51.100.0/24 1; 2001:db8::/32 1;")

  out3 <- withr::local_tempdir()
  withr::with_envvar(c(ADMIN_ALLOW_CIDR = "192.0.2.10"), expect_identical(render(out3), 0L))
  expect_identical(geo_line(out3), "192.0.2.10 1;")

  # Anything but addresses would be nginx syntax written into the config.
  out4 <- withr::local_tempdir()
  for (bad in c("1.2.3.4;deny", "all", "1.2.3.4/24}", "example.org")) {
    expect_false(identical(render(out4, "--admin-allow", shQuote(bad)), 0L), info = bad)
  }
  expect_length(list.files(out4, recursive = TRUE), 0L)
  expect_true(any(grepl("^ADMIN_ALLOW_CIDR=", read_deploy(root, "host.env.template"))))
})

test_that("the services nginx fronts listen on loopback only", {
  root <- skip_unless_deploy()
  compose <- read_deploy(root, "keycloak", "docker-compose.yml.template")
  published <- grep('^\\s*- "[0-9.:]+:[0-9]+"', compose, value = TRUE)
  expect_gte(length(published), 1L)
  expect_true(all(grepl('- "127\\.0\\.0\\.1:', published)))
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  expect_identical(yaml_value(sp, "target-bind-ip"), "127.0.0.1")
})

# ---- container confinement ---------------------------------------------------

test_that("app containers are confined with the settings ShinyProxy supports", {
  root <- skip_unless_deploy()
  sp <- uncommented(read_deploy(root, "shinyproxy", "application.yml.template"))
  expect_identical(yaml_value(sp, "container-privileged"), "false")
  expect_match(yaml_value(sp, "container-memory-limit"), "^[0-9]+[mg]$")
  expect_match(yaml_value(sp, "container-cpu-limit"), "^[0-9.]+$")
  expect_identical(yaml_value(sp, "container-network"), "omicsapp-net")
  # Keys ShinyProxy's Docker backend does not have would be inert --
  # hardening on paper only. They belong elsewhere (README).
  for (key in c("container-cap-drop", "cap-drop", "container-read-only", "read-only",
                "container-pids-limit", "pids-limit", "security-opt",
                "container-security-opt", "no-new-privileges")) {
    expect_true(is.na(yaml_value(sp, key)), info = key)
  }
})

test_that("the image gives its user no way up and caps its processes", {
  root <- skip_unless_deploy()
  docker <- uncommented(read_deploy(root, "docker", "Dockerfile"))
  expect_true(any(grepl("-perm -4000 -o -perm -2000 \\) -exec chmod ug-s", docker, fixed = TRUE)))
  expect_true(any(grepl("setuid/setgid files remain", docker, fixed = TRUE)))
  users <- sub("^USER\\s+", "", grep("^USER", docker, value = TRUE))
  expect_identical(utils::tail(users, 1L), "omics")
  cmd <- grep("^CMD", docker, value = TRUE)
  expect_match(cmd, "ulimit -S -u", fixed = TRUE)
  expect_match(cmd, "OMICSAPP_MAX_PROCS", fixed = TRUE)
  expect_match(cmd, "exec R", fixed = TRUE)
})

test_that("the compose services are confined, rotated and limited", {
  root <- skip_unless_deploy()
  compose <- read_deploy(root, "keycloak", "docker-compose.yml.template")
  last <- grep("^networks:", compose) - 1L
  starts <- grep("^  [a-z-]+:\\s*$", compose)
  starts <- starts[starts > grep("^services:", compose) & starts < last]
  ends <- c(starts[-1L] - 1L, last)
  services <- Map(function(s, e) compose[s:e], starts, ends)
  names(services) <- trimws(sub(":.*", "", compose[starts]))
  expect_setequal(names(services), c("keycloak-db", "keycloak"))
  for (nm in names(services)) {
    svc <- uncommented(services[[nm]])
    expect_true(any(grepl('security_opt: \\["no-new-privileges:true"\\]', svc)), info = nm)
    expect_true(any(grepl("cap_drop: \\[ALL\\]", svc)), info = nm)
    expect_true(any(grepl("pids_limit: [0-9]+", svc)), info = nm)
    expect_true(any(grepl("mem_limit: [0-9]+[mg]", svc)), info = nm)
    expect_true(any(grepl("cpus: [0-9.]+", svc)), info = nm)
    expect_true(any(grepl('max-size: "[0-9]+m"', svc)), info = nm)
  }
  db <- uncommented(services[["keycloak-db"]])
  expect_true(any(grepl("read_only: true", db, fixed = TRUE)))
  expect_true(any(grepl("- /var/run/postgresql", db, fixed = TRUE)))
})

test_that("the daemon defaults are valid JSON with rotation and no-new-privileges", {
  root <- skip_unless_deploy()
  skip_if_not_installed("jsonlite")
  d <- jsonlite::fromJSON(file.path(root, "docker", "daemon.json"))
  expect_true(isTRUE(d[["no-new-privileges"]]))
  expect_false(is.null(d[["log-opts"]][["max-size"]]))
})

# ---- egress -------------------------------------------------------------------

test_that("the app network is the one egress.sh and the README name", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  net <- yaml_value(sp, "container-network")
  egress <- read_deploy(root, "scripts", "egress.sh")
  expect_true(any(grepl(sprintf('NETWORK="${NETWORK:-%s}"', net), egress, fixed = TRUE)))
  readme <- read_deploy(root, "README.md")
  expect_true(any(grepl(sprintf("docker network create -o com.docker.network.bridge.name=[a-z-]+ %s$", net),
                        readme)))
  expect_true(any(grepl("deploy/scripts/egress.sh apply", readme, fixed = TRUE)))
  unit <- read_deploy(root, "systemd", "omicsapp-egress.service")
  expect_true(any(grepl("^ExecStart=.*/deploy/scripts/egress.sh apply$", unit)))
  expect_true(any(grepl("^After=docker.service$", unit)))
})

test_that("egress.sh drops new connections from the app bridge, idempotently", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  t <- withr::local_tempdir()
  state <- file.path(t, "rules")
  file.create(state)
  writeLines(c("#!/bin/sh",
               "case \"$*\" in *'{{index .Options'*) echo br-omicsapp ;; *) echo deadbeef0000ffff ;; esac"),
             file.path(t, "docker"))
  # A fake iptables keeping its rules in a file: -n -L, -C, -I, -D.
  writeLines(c("#!/bin/bash",
               sprintf("STATE=%s", shQuote(state)),
               "op=\"$1\"; shift",
               "case \"$op\" in",
               "  -n) exit 0 ;;",
               "  -C) grep -qxF -- \"$*\" \"$STATE\" ;;",
               "  -I) chain=\"$1\"; shift 2; echo \"$chain $*\" >> \"$STATE\" ;;",
               "  -D) grep -vxF -- \"$*\" \"$STATE\" > \"$STATE.t\"; mv \"$STATE.t\" \"$STATE\" ;;",
               "esac"), file.path(t, "iptables"))
  Sys.chmod(file.path(t, c("docker", "iptables")), "755")
  run <- function(action) {
    withr::with_envvar(c(DOCKER = file.path(t, "docker"), IPTABLES = file.path(t, "iptables"),
                         IP6TABLES = file.path(t, "no-ip6tables")),
      system2("bash", c(file.path(root, "scripts", "egress.sh"), action), stdout = FALSE, stderr = FALSE))
  }
  expect_false(run("check") == 0L)
  expect_identical(run("apply"), 0L)
  expect_identical(run("apply"), 0L)
  rules <- readLines(state)
  expect_length(rules, 2L)
  expect_true(any(grepl("^DOCKER-USER -i br-omicsapp -m conntrack --ctstate NEW,INVALID .* -j DROP$", rules)))
  expect_true(any(grepl("^INPUT -i br-omicsapp -m conntrack --ctstate NEW,INVALID .* -j DROP$", rules)))
  expect_identical(run("check"), 0L)
  expect_identical(run("remove"), 0L)
  expect_length(readLines(state), 0L)
})

test_that("app containers do not try to refresh gene sets over the network", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  expect_identical(yaml_value(sp, "OMICSCORE_GENESET_TTL_DAYS"), '"0"')
})

# ---- logs: kept on the host, rotated, backed up --------------------------------

test_that("every log is on the host, rotated, and in the backup", {
  root <- skip_unless_deploy()
  sp <- read_deploy(root, "shinyproxy", "application.yml.template")
  container_logs <- yaml_value(sp, "container-log-path")
  expect_match(container_logs, "^/var/log/shinyproxy/")
  expect_match(yaml_value(sp, "max-file-size"), "^[0-9]+MB$")
  expect_match(yaml_value(sp, "max-history"), "^[0-9]+$")
  expect_match(yaml_value(sp, "total-size-cap"), "^[0-9]+[MG]B$")

  rotate <- read_deploy(root, "logrotate", "omicsapp")
  expect_true(any(grepl(paste0("^", container_logs, "/\\*\\.log \\{"), rotate)))
  # ShinyProxy rotates its own file; logrotate must not do it too.
  expect_false(any(grepl("shinyproxy.log", uncommented(rotate), fixed = TRUE)))

  nginx <- uncommented(read_deploy(root, "nginx", "omicsapp.conf.template"))
  expect_true(any(grepl("access_log /var/log/nginx/omicsapp\\.access\\.log;", nginx)))

  backup <- read_deploy(root, "scripts", "backup.sh")
  log_paths <- grep('^LOG_PATHS=', backup, value = TRUE)
  expect_length(log_paths, 1L)
  for (p in c("/var/log/shinyproxy", "/var/log/nginx", "/var/log/omicsapp-backup.log")) {
    expect_match(log_paths, p, fixed = TRUE)
  }
  extra <- grep('^EXTRA_PATHS=', backup, value = TRUE)
  for (p in c("/etc/docker/daemon.json", "/etc/logrotate.d/omicsapp",
              "/etc/systemd/system/omicsapp-egress.service")) {
    expect_match(extra, p, fixed = TRUE)
  }
  expect_true(any(grepl("^LOG_PATHS=", read_deploy(root, "backup.env.template"))))
})

# ---- rollback: immutable tags, a script to switch, digests ---------------------

fake_bin <- function(dir, name, lines) {
  path <- file.path(dir, name)
  writeLines(c("#!/bin/bash", lines), path)
  Sys.chmod(path, "755")
  path
}

test_that("build_image.sh tags version and commit, and refuses latest", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  t <- withr::local_tempdir()
  calls <- file.path(t, "calls")
  docker <- fake_bin(t, "docker", c(
    sprintf("echo \"$*\" >> %s", shQuote(calls)),
    "case \"$1 $2\" in 'image inspect') exit 1 ;; esac", "exit 0"))
  build <- file.path(root, "scripts", "build_image.sh")
  status <- withr::with_envvar(c(DOCKER = docker),
    system2("bash", c(build, "omicsapp:1.0"), stdout = FALSE, stderr = FALSE))
  expect_identical(status, 0L)
  version <- read.dcf(file.path(root, "..", "packages", "omicsApp", "DESCRIPTION"))[1, "Version"]
  built <- grep("^build ", readLines(calls), value = TRUE)
  expect_length(built, 1L)
  tags <- regmatches(built, gregexpr("-t [^ ]+", built))[[1L]]
  immutable <- sub("^-t ", "", tags[[1L]])
  expect_match(immutable, sprintf("^omicsapp:%s-([0-9a-f]{12}(-dirty\\.[0-9]{14})?|norev\\.[0-9]{14})$",
                                  gsub(".", "\\.", version, fixed = TRUE)))
  expect_true("-t omicsapp:1.0" %in% tags)
  expect_match(built, "org.opencontainers.image.revision=", fixed = TRUE)

  # `latest`, alone or as an alias, is refused and nothing is built.
  file.remove(calls)
  for (bad in c("omicsapp:latest", "latest")) {
    status <- withr::with_envvar(c(DOCKER = docker),
      system2("bash", c(build, bad), stdout = FALSE, stderr = FALSE))
    expect_false(identical(status, 0L), info = bad)
  }
  expect_false(file.exists(calls) && any(grepl("^build ", readLines(calls))))

  # An immutable tag that already exists is never rebuilt over.
  docker2 <- fake_bin(t, "docker2", c(sprintf("echo \"$*\" >> %s", shQuote(calls)), "exit 0"))
  withr::with_envvar(c(DOCKER = docker2),
    system2("bash", build, stdout = FALSE, stderr = FALSE))
  expect_false(any(grepl("^build ", readLines(calls))))
})

test_that("rollback.sh switches the image, keeps the old config, and restarts", {
  root <- skip_unless_deploy()
  skip_unless_bash()
  t <- withr::local_tempdir()
  conf <- file.path(t, "application.yml")
  file.copy(file.path(root, "shinyproxy", "application.yml.template"), conf)
  Sys.chmod(conf, "600")
  restarts <- file.path(t, "restarts")
  docker <- fake_bin(t, "docker",
    "if [ \"$1 $2\" = 'image inspect' ]; then case \"$3\" in omicsapp:old-*|omicsapp:new-*) exit 0 ;; *) exit 1 ;; esac; fi")
  systemctl <- fake_bin(t, "systemctl", sprintf("echo \"$*\" >> %s", shQuote(restarts)))
  fake_bin(t, "curl", "exit 0")
  script <- file.path(root, "scripts", "rollback.sh")
  run <- function(...) {
    withr::with_envvar(c(SP_CONF = conf, DOCKER = docker, SYSTEMCTL = systemctl,
                         WAIT_SECONDS = "2", PATH = paste(t, Sys.getenv("PATH"), sep = ":")),
      system2("bash", c(script, ...), stdout = FALSE, stderr = FALSE, stdin = "/dev/null"))
  }
  image <- function() yaml_value(readLines(conf), "container-image")
  before <- image()

  expect_identical(run("--yes", "omicsapp:new-1a2b3c4d5e6f"), 0L)
  expect_identical(image(), "omicsapp:new-1a2b3c4d5e6f")
  expect_identical(readLines(restarts), "restart shinyproxy")
  history <- readLines(file.path(t, "image-history"))
  expect_match(history, sprintf("%s -> omicsapp:new-1a2b3c4d5e6f$", before))
  backups <- list.files(t, pattern = "^application\\.yml\\.bak-")
  expect_length(backups, 1L)
  expect_identical(yaml_value(readLines(file.path(t, backups)), "container-image"), before)
  # Only that one line changed; comments and everything else are intact.
  changed <- readLines(conf) != readLines(file.path(t, backups))
  expect_identical(sum(changed), 1L)
  if (.Platform$OS.type == "unix") expect_identical(format(file.info(conf)$mode), "600")

  # And back.
  Sys.sleep(1.1)   # a distinct .bak name
  expect_identical(run("--yes", "omicsapp:old-0f0f0f0f0f0f"), 0L)
  expect_identical(image(), "omicsapp:old-0f0f0f0f0f0f")
  expect_length(readLines(file.path(t, "image-history")), 2L)

  # Refused, and nothing touched: an image not on the host, `latest`, no
  # confirmation without a terminal.
  for (args in list(c("--yes", "omicsapp:missing-123"), c("--yes", "omicsapp:latest"),
                    "omicsapp:new-1a2b3c4d5e6f")) {
    expect_false(identical(run(args), 0L), info = paste(args, collapse = " "))
    expect_identical(image(), "omicsapp:old-0f0f0f0f0f0f")
  }
  expect_length(readLines(restarts), 2L)
})

test_that("backup and the restore drill follow the image ShinyProxy runs", {
  root <- skip_unless_deploy()
  for (script in c("backup.sh", "restore_check.sh")) {
    sh <- read_deploy(root, "scripts", script)
    expect_true(any(grepl("container-image", sh, fixed = TRUE)), info = script)
    expect_true(any(grepl('SP_CONF="${SP_CONF:-/etc/shinyproxy/application.yml}"', sh, fixed = TRUE)),
                info = script)
  }
  expect_true("APP_IMAGE=" %in% read_deploy(root, "backup.env.template"))
})

test_that("CI tags the image with version and commit, and smoke-tests it confined", {
  root <- skip_unless_deploy()
  wf <- file.path(dirname(root), ".github", "workflows", "production-image.yaml")
  skip_if(!file.exists(wf), "workflows not present")
  lines <- readLines(wf, warn = FALSE)
  expect_true(any(grepl("immutable=omicsapp:${version}-${GITHUB_SHA:0:12}", lines, fixed = TRUE)))
  expect_true(any(grepl("${{ steps.tag.outputs.immutable }}", lines, fixed = TRUE)))
  expect_true(any(grepl("org.opencontainers.image.revision=${{ github.sha }}", lines, fixed = TRUE)))
  expect_true(any(grepl("--read-only", lines, fixed = TRUE)))
  expect_true(any(grepl("--cap-drop ALL --security-opt no-new-privileges", lines, fixed = TRUE)))
})

# Every third-party image reference, as written.
image_references <- function(root) {
  docker <- read_deploy(root, "docker", "Dockerfile")
  compose <- read_deploy(root, "keycloak", "docker-compose.yml.template")
  c(sub("^ARG BASE_IMAGE=", "", grep("^ARG BASE_IMAGE=", docker, value = TRUE)),
    trimws(sub("^\\s+image:\\s*", "", grep("^\\s+image:", compose, value = TRUE))))
}

test_that("the Dockerfile builds FROM the pinnable BASE_IMAGE argument", {
  root <- skip_unless_deploy()
  docker <- uncommented(read_deploy(root, "docker", "Dockerfile"))
  from <- grep("^FROM ", docker)
  arg <- grep("^ARG BASE_IMAGE=", docker)
  expect_length(from, 1L)
  expect_length(arg, 1L)
  expect_lt(arg, from)
  expect_identical(docker[[from]], "FROM ${BASE_IMAGE}")
})

test_that("every image with a recorded digest is referenced by that digest", {
  root <- skip_unless_deploy()
  lock_file <- file.path(root, "docker", "base-digests.lock")
  lock <- if (file.exists(lock_file)) uncommented(readLines(lock_file, warn = FALSE)) else character(0)
  lock <- lock[nzchar(trimws(lock))]
  recorded <- stats::setNames(sub("^\\S+\\s+", "", lock), sub("\\s.*$", "", lock))
  expect_true(all(grepl("^sha256:[0-9a-f]{64}$", recorded)))
  refs <- image_references(root)
  expect_gte(length(refs), 3L)
  unpinned <- character(0)
  for (ref in refs) {
    name <- sub("@sha256:.*$", "", ref)
    if (name %in% names(recorded)) {
      # Once a digest is recorded, a reference without it -- or with a
      # different one -- is a build that no longer matches the record.
      expect_identical(ref, paste0(name, "@", recorded[[name]]), info = name)
    } else {
      unpinned <- c(unpinned, name)
    }
  }
  # Not yet a failure: a digest can only be recorded with access to the
  # registry (see the TODO in deploy/README.md).
  if (length(unpinned)) {
    warning("Pinned by tag only, no digest recorded: ", paste(unpinned, collapse = ", "),
            ". Run deploy/scripts/pin_base_digests.sh with registry access.", call. = FALSE)
  }
  # Nothing recorded that is no longer used.
  expect_length(setdiff(names(recorded), sub("@sha256:.*$", "", refs)), 0L)
})
