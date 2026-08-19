# Deploying omicsApp

Serves the app to a small group (3–4 concurrent users) on an internal
network, one container per logged-in user.

```
browser ──HTTP──> nginx ──> ShinyProxy(:8080) ──> omicsapp container × N
                              │ simple auth        └─ /data ← /srv/omicsapp/users/<user>
                              └─ lifecycle
```

Per-user containers are the point: each user's expression matrices and
results live only in their own process, and one session running out of
memory cannot take anyone else down with it.

## Contents

| Path | What it is |
|---|---|
| `docker/Dockerfile` | The image. Build from the **repository root**. |
| `docker/prewarm_genesets.R` | Bakes the MSigDB tables in at build time. |
| `shinyproxy/application.yml.template` | Copy to `application.yml` on the server and fill in users. |
| `nginx/omicsapp.conf` | Reverse proxy, including the WebSocket headers Shiny needs. |
| `scripts/build_image.sh` | Wrapper that gets the build context right. |
| `scripts/add_user.sh` | Creates a user's private directory with the ownership the container needs. |

## First deployment

Verify each step before starting the next. Standing the whole stack up
and then finding it broken leaves five candidate causes; this order
localises a failure to the step that caused it, and puts the two
riskiest things first so an hour is not spent before finding out.

### 0. Decide what you are deploying

The Dockerfile COPYs the **working tree**, not a commit, so whatever is
checked out is what ships.

```bash
git checkout main
git status --short   # must be empty, or the image matches no commit
```

### 1. Pre-flight (2 minutes, saves an hour)

The base image tag is pinned in the Dockerfile but is worth confirming
against what your host can actually pull, along with the R and
Bioconductor versions it carries:

```bash
docker run --rm bioconductor/bioconductor_docker:RELEASE_3_20 \
  R -q -e 'cat(R.version.string, "| Bioc", as.character(BiocManager::version()), "\n")'
df -h /var/lib/docker   # image is 5-7 GB; build cache wants 2-3x that
```

### 2. Build (30-60 minutes)

```bash
deploy/scripts/build_image.sh omicsapp:1.0
```

A failure part way through is not a restart: Docker caches each layer,
so a fix re-runs only from the layer that failed. The Bioconductor and
LaTeX layers are the likely ones. Dropping the LaTeX layer costs the
PDF report and saves 1.5 GB.

### 3. Smoke-test the container on its own — do not skip

This separates "does the app work" from "is ShinyProxy configured
right". Debugging both at once is what makes a deployment take a day.

```bash
docker run --rm -p 3838:3838 -e OMICSAPP_DATA_DIR=/tmp/data omicsapp:1.0
```

Open `http://<host>:3838` and confirm all seven views render and the
Report view shows an "Analysis code" card. Then check the two things
that fail silently rather than loudly:

```bash
# Gene-set cache is in the image. Missing, the app still works -- it
# just pays ~10s on every user's first enrichment, forever.
docker run --rm omicsapp:1.0 \
  R -q -e 'cat(length(list.files(Sys.getenv("OMICSCORE_GENESET_CACHE"))), "cached tables\n")'
# expect 14

# The future plan is parallel. Sequential means one user's DESeq2
# freezes every other session, with nothing on screen to say so.
docker run --rm omicsapp:1.0 \
  R -q -e 'library(future); cat(class(future::plan())[2], "\n")'
```

### 4. Network and storage

```bash
docker network create sp-net
sudo mkdir -p /srv/omicsapp/users
sudo deploy/scripts/add_user.sh puweilin   # once per user
```

### 5. ShinyProxy

```bash
sudo mkdir -p /etc/shinyproxy
sudo cp deploy/shinyproxy/application.yml.template /etc/shinyproxy/application.yml
sudo chmod 600 /etc/shinyproxy/application.yml
htpasswd -bnBC 10 "" 'the-password' | tr -d ':\n'   # paste into proxy.users
```

Log in with one test account before hashing everyone's password: the
`{bcrypt}` prefix depends on the ShinyProxy version.

### 6. nginx and firewall

```bash
sudo cp deploy/nginx/omicsapp.conf /etc/nginx/sites-available/omicsapp
sudo ln -s /etc/nginx/sites-available/omicsapp /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx
sudo ufw allow from 192.168.51.0/24 to any port 80 proto tcp
```

Users then reach the app at `http://192.168.51.52`.

### 7. Acceptance

Walk one real dataset through, checking each line:

| Action | Expected |
|---|---|
| Upload, confirm import | Project appears |
| Re-upload the **same** file | "already loaded" — nothing is cleared |
| Upload a **different** file | Dialog naming the analyses that will be cleared |
| Run QC and a differential analysis | Volcano is two-coloured; caption reads `adj_p_value < 0.05` |
| Drag the FDR slider | Hit table changes, **volcano does not** |
| Report view | "Analysis code" shows the calls; downloads as `.R` |
| Project view, Save as | Appears under "My projects" |
| Close the tab, log back in | "Restore last session" works |
| `ls /srv/omicsapp/users/<user>/` | `.omp` files and a `raw/` directory |
| **Two people run an analysis at once** | Neither waits for the other |

The last row is the Phase 0 fix that mattered most; it needs a
colleague to check.

## Things that will bite you

**The build context is the repository root.** The Dockerfile COPYs both
`packages/omicsCore` and `packages/omicsApp`. `docker build` from
`deploy/docker/` fails on those lines. Use `build_image.sh`, or pass
`-f deploy/docker/Dockerfile .` yourself.

**Create a user's directory before their first login.** Docker creates a
missing bind-mount source owned by `root`; the container runs as uid
1001 and then cannot write, so the user silently loses the ability to
save anything. `add_user.sh` handles this.

**Never commit `application.yml`.** `.gitignore` excludes it. Git history
is permanent, and a private repository still gets cloned to laptops.

**ShinyProxy has no idle detection.** A tab left open counts as active
indefinitely. `heartbeat-timeout` only catches a *closed* tab;
`default-proxy-max-lifetime` is the only thing that reclaims a forgotten
one.

**Restarting ShinyProxy drops running sessions.** Adding a user requires
a restart, so do it when nobody is mid-analysis.

**Check the ShinyProxy version against the template.**
`minimum-seats-available` is 3.x syntax (2.x used
`container-pre-initialization`), and the `{bcrypt}` password prefix is
worth confirming with one test login before hashing everyone's password.

## Resources

| | |
|---|---|
| RAM | ≥ 32 GB (4 × 6 GB container limit, plus gene sets resident per container) |
| CPU | ≥ 8 cores; `OMP_NUM_THREADS=1` in the image stops workers oversubscribing |
| Disk | 1 TB SSD is comfortable. Growth is ~50–75 GB/year at three analyses a day |
| GPU | not used |
| Image | 5–7 GB (≈1.5 GB of that is the LaTeX layer for PDF reports) |

## Environment variables the image honours

| Variable | Default | Meaning |
|---|---|---|
| `OMICSAPP_DATA_DIR` | `/data` | Per-user project directory. Falls back to a temp dir when unset, so nothing writes to a home directory by surprise. |
| `OMICSAPP_QUOTA_GB` | `50` | Soft quota. Invalid or unset means unlimited — a typo in the config must not lock a user out of their own data. |
| `OMICSCORE_GENESET_CACHE` | `/opt/genesets` | Pre-built MSigDB tables. A missing or unreadable cache costs ~10s on first enrichment, never a wrong answer. |
| `OMP_NUM_THREADS` | `1` | One BLAS thread per worker. |

## Operations

**Backup.** `/srv/omicsapp/users` is the whole of it. `.omp` files are
the irreplaceable part — raw uploads are usually still obtainable from
the instrument, but the analysis parameters encoded in a project are
not.

```bash
rsync -a --delete /srv/omicsapp/users/ /backup/omicsapp-users/
```

**Rebuild after a code change.** Layers up to the package COPY are
cached, so a rebuild is a couple of minutes rather than an hour.

```bash
deploy/scripts/build_image.sh omicsapp:1.1
# then update container-image in application.yml and restart ShinyProxy
```

**Where the logs are.** ShinyProxy writes to
`/var/log/shinyproxy/shinyproxy.log`; a container's own R output is in
`docker logs <container>`.
