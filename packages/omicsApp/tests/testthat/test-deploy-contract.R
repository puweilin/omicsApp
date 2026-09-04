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
  nginx <- read_deploy(root, "nginx", "omicsapp.conf")
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
  nginx <- read_deploy(root, "nginx", "omicsapp.conf")
  code <- nginx[!grepl("^\\s*#", nginx)]
  expect_false(any(grepl("http2", code)))
  expect_true(any(grepl("proxy_set_header Upgrade\\s+\\$http_upgrade", code)))
  expect_true(any(grepl("proxy_set_header Connection\\s+\"upgrade\"", code)))
  expect_gte(sum(grepl("X-Forwarded-Proto\\s+\\$scheme", code)), 2L)
  expect_true(any(grepl("proxy_read_timeout", code)))
})

test_that("nginx accepts uploads at least as large as the app does", {
  root <- skip_unless_deploy()
  nginx <- read_deploy(root, "nginx", "omicsapp.conf")
  line <- grep("client_max_body_size", nginx, value = TRUE)
  expect_length(line, 1L)
  nginx_mb <- as.numeric(sub(".*client_max_body_size\\s+([0-9]+)M;.*", "\\1", line))
  app_mb <- formals(launch)$max_upload_mb
  expect_gte(nginx_mb, app_mb)
})

# ---- the password floor, which drifted across three files -----------------

test_that("the password policy is the same number everywhere it is written", {
  root <- skip_unless_deploy()
  realm <- read_deploy(root, "keycloak", "omicsapp-realm.json")
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
