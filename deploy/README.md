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
sudo mkdir -p /srv/omicsapp/users            # SSD
sudo mkdir -p /mnt/hdd/omicsapp/raw          # HDD, if splitting
sudo OMICSAPP_RAW_ROOT=/mnt/hdd/omicsapp/raw \
     deploy/scripts/add_user.sh puweilin     # once per user
```

See [Storage layout](#storage-layout) for what belongs on which disk.
Leave `OMICSAPP_RAW_ROOT` unset to keep everything on one.

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

**Create a user's directories before their first login.** Docker creates
a missing bind-mount source owned by `root`; the container runs
unprivileged and then cannot write, so the user silently loses the
ability to save anything — the app reports a failed save and carries on.
`add_user.sh` handles this.

**The container uid and the directory owner must match.** The kernel
compares numbers, not names: `omics` inside the image and any account on
the host are unrelated. The default is 1001 on both sides. If you built
with `APP_UID=$(id -u)` because you have no root, pass the same value to
`add_user.sh`, or the two halves will disagree and produce exactly the
silent failure above.

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

## Storage layout

Measured, on a 20k feature by 60 sample workbook: a `.omp` project file
is **5.9 MB**, the raw upload it came from is **20 MB**, and one full
analysis produces about **45 MB** all told. At three analyses a day that
is **50–75 GB a year**.

| What | Growth | Access | Disk |
|---|---|---|---|
| `.omp` project files | ~35 GB/yr | Read and written on every save, open, and finished analysis | **SSD** |
| `raw/` archived uploads | ~15 GB/yr | Written once, almost never read | **HDD** |
| Docker images and layers | 7 GB, fixed | Read on every container start | **SSD, not negotiable** |
| Backups | mirrors the above | Nightly | **HDD** |

Two things worth knowing before moving anything:

**There is no capacity pressure on the SSD.** 75 GB a year against 4 TB
is fifty years. Splitting the archive onto the HDD is about putting bulk
where bulk belongs, not about running out of room.

**The latency cost of the HDD is small but not zero.** A `.omp` write is
0.06 s on SSD and most of that is serialisation rather than I/O; a
spinning disk adds roughly 40–70 ms. Autosave fires five to seven times
per workflow, so projects on the HDD would be felt slightly. Archived
uploads are written once and never re-read, so they cost nothing there.

**`/var/lib/docker` must stay on the SSD.** The image is 7 GB and every
container start reads from it; a build is heavy random I/O. This is the
one placement that is not a preference.

One clarification on the 7 GB, because it changes the arithmetic: image
layers are **read-only and shared**. Four concurrent containers do not
use 4 × 7 GB — there is one copy on disk, and each container adds only
its writable layer, which stays tiny here because all data goes to bind
mounts rather than into the container.

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
| `OMICSCORE_GENESET_CACHE` | `/opt/genesets` | Pre-built gene-set tables (MSigDB at build time, current KEGG after a refresh — see Operations). A missing or unreadable cache costs ~10s on first enrichment, never a wrong answer. |
| `OMICSCORE_GENESET_TTL_DAYS` | `30` | Age at which a **live-sourced** KEGG cache re-fetches itself from KEGG REST on next use. `0` disables. Prewarmed MSigDB tables never trigger network calls. |
| `OMP_NUM_THREADS` | `1` | One BLAS thread per worker. |

## Operations

**Backup — do this on day one.** Losing a disk is far likelier than any
attack, and `.omp` files are the irreplaceable part: a raw upload can
usually be exported from the instrument again, the analysis parameters
encoded in a project cannot.

```bash
# Nightly, e.g. from cron. The 60 TB HDD is the obvious target.
rsync -a --delete /srv/omicsapp/users/ /mnt/hdd/omicsapp-backup/users/
rsync -a --delete /mnt/hdd/omicsapp/raw/ /mnt/hdd/omicsapp-backup/raw/
```

Backing the archive up onto the same disk it lives on protects against
deletion but not against that disk failing. If the raw archive matters
to you, send it somewhere else as well.

**Keeping KEGG current.** The image bakes the MSigDB tables in, and the
`kegg` table among them is the 2011 `KEGG_LEGACY` snapshot (186 human
pathways; KEGG today has ~370). `omicsCore::refresh_geneset_cache()`
replaces it with a live fetch from the KEGG REST API and re-snapshots
the other databases from msigdbr. To keep that current without
rebuilding images, move the cache onto a volume and refresh it on a
schedule:

```bash
# one-time: seed the volume, then mount it in application.yml with
#   container-volumes: [ ..., "/srv/omicsapp/genesets:/opt/genesets" ]
mkdir -p /srv/omicsapp/genesets
docker run --rm -v /srv/omicsapp/genesets:/opt/genesets omicsapp:1.0 \
  R -q -e 'for (org in c("Hs","Mm")) omicsCore::refresh_geneset_cache(organism = org, force = TRUE)'
```

```bash
# /etc/cron.d/omicsapp-genesets — monthly, 03:00 on the 1st
0 3 1 * * root docker run --rm -v /srv/omicsapp/genesets:/opt/genesets omicsapp:1.0 R -q -e 'for (org in c("Hs","Mm")) omicsCore::refresh_geneset_cache(organism = org, force = TRUE)' >> /var/log/omicsapp-genesets.log 2>&1
```

The mount shadows the baked copy, which is why the seed run writes all
databases, not just KEGG. A failed fetch keeps the previous file, so
the worst case of a dead network on refresh night is a month-old table.
Every result records which definitions it used
(`bundle$params$geneset_sources`), and `omicsCore::geneset_cache_status()`
shows what is live on disk. One caution: KEGG's license permits this
per-query use but not redistribution — treat the volume as
deployment-local data, never as something to publish or commit.

**Rebuild after a code change.** Layers up to the package COPY are
cached, so a rebuild is a couple of minutes rather than an hour.

```bash
deploy/scripts/build_image.sh omicsapp:1.1
# then update container-image in application.yml and restart ShinyProxy
```

**Where the logs are.** ShinyProxy writes to
`/var/log/shinyproxy/shinyproxy.log`; a container's own R output is in
`docker logs <container>`.

## Hardening, once it works

Leave these until the acceptance checklist passes. First deployments go
wrong for ordinary reasons, and every extra variable is one more
candidate.

**Give ShinyProxy its own account.** The official package runs it as
root. It does not need to be — it needs to reach the Docker socket,
which is group membership, not privilege:

```bash
sudo useradd --system --no-create-home --shell /usr/sbin/nologin shinyproxy
sudo usermod -aG docker shinyproxy
sudo chown shinyproxy /etc/shinyproxy/application.yml
# then, in the systemd unit:  User=shinyproxy   Group=docker
```

A dedicated account rather than a person's: a stolen SSH key should not
hand over the service, and a compromised service should not reach
someone's files.

Be clear about what this buys, though. **Membership of the `docker`
group is equivalent to root** — anything that can reach the socket can
ask the daemon to start a privileged container mounting the host's
filesystem. So this raises the cost of a ShinyProxy vulnerability by a
step; it is not a boundary.

The boundary is one layer down, and it is already in place: containers
run unprivileged, and the Docker socket is not mounted into any of them.
Those two are what a compromise in one of the image's 224 R packages
would run into. Keep them, and rebuild the image periodically so those
packages get their patches — dependencies age whether or not anyone is
watching.
