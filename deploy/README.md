# Deploying omicsApp

Serves the app to a small group (3–4 concurrent users) on an internal
network, one container per logged-in user.

```
browser ──HTTPS──> nginx ─┬─> ShinyProxy(:8081) ──> omicsapp container × N
                          │     └─ lifecycle          └─ /data ← /srv/omicsapp/users/<sub>
                          └─> Keycloak(:8180) ──> accounts, passwords, self-service
```

ShinyProxy holds no credentials. Keycloak owns them, which is what lets
a user change their own password and username without an admin editing a
file and restarting a service that would drop everyone's running
analysis. Storage is keyed on the Keycloak `sub` — an immutable UUID —
so a rename does not orphan anyone's projects. See
[`keycloak/README.md`](keycloak/README.md).

Per-user containers are the point: each user's expression matrices and
results live only in their own process, and one session running out of
memory cannot take anyone else down with it.

## Contents

| Path | What it is |
|---|---|
| `docker/Dockerfile` | The image. Build from the **repository root**. |
| `docker/prewarm_genesets.R` | Bakes the MSigDB tables in at build time. |
| `shinyproxy/application.yml.template` | Copy to `application.yml` on the server and fill in the client secret. |
| `nginx/omicsapp.conf` | Reverse proxy: TLS, the WebSocket headers Shiny needs, and Keycloak at `/auth`. |
| `keycloak/` | The identity provider: compose file, realm definition, and its own README. |
| `scripts/build_image.sh` | Wrapper that gets the build context right. |
| `scripts/add_user.sh` | Creates an account in Keycloak and its storage directory, in that order. |
| `scripts/list_users.sh` | Maps the UUID directory names back to people. |

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

**Getting it onto the server.** There are no git credentials there, so
this is a copy from the machine you work on, not a `git pull`:

```bash
rsync -avz --delete \
  --exclude '.git/' \
  --exclude 'deploy/keycloak/.env' \
  ~/SciProject/CHISSS/omicsApp/ ps@192.168.51.52:~/omicsApp/
```

Both flags earn their place. `--delete` is what removes a script that
was deleted upstream — otherwise it lingers on the server and someone
runs it a year later. The exclusion protects `deploy/keycloak/.env`,
which holds the Keycloak secrets, exists only on the server, and would
otherwise be deleted by the next sync. (rsync does not delete excluded
files, so naming it here is what keeps it.)

### 1. Host prerequisites

Four things, and only four — everything else the application needs is
inside the image.

```bash
sudo apt-get update
sudo apt-get install -y docker.io nginx rsync
sudo systemctl enable --now docker

# ShinyProxy is a jar; the .deb installs it plus a systemd unit.
#
# 3.2.4, not 3.1.1: on Docker 28+ the older release cannot start any
# container at all. See "Things that will bite you".
#
# openjdk-21-jdk-headless, not -jre-headless: the .deb depends on
# "openjdk-21-jdk-headless | openjdk-21-jre", and -jre-headless is
# neither of those -- it is the package they both depend on. dpkg
# unpacks and then refuses to configure, which reads like a broken
# download rather than a missing dependency. (3.1.1 wanted Java 17;
# 3.2.x wants 21.)
sudo apt-get install -y openjdk-21-jdk-headless
wget https://github.com/openanalytics/shinyproxy/releases/download/v3.2.4/shinyproxy_3.2.4_amd64.deb
sudo dpkg -i shinyproxy_3.2.4_amd64.deb

# Optional: the mDNS hostname alias (see nginx/omicsapp.conf)
sudo apt-get install -y avahi-daemon
sudo hostnamectl set-hostname omics
```

Confirm before continuing:

```bash
docker run --rm hello-world     # daemon is up and you can reach it
docker version --format '{{.Server.Version}}'
systemctl status shinyproxy     # installed; will fail to start until configured
nginx -v
```

**On a host older than 24.04, look at that engine version.** The
application image is built on Ubuntu 24.04 (see below), so on an older
host its glibc is newer than the host's — fine in itself, except that
Docker engines before 20.10.10 shipped a seccomp profile that blocked
`clone3`, a syscall glibc 2.34+ uses. The symptom is not a clear error
but a container that dies on start. If `docker.io` gives you anything
older, install from Docker's own repository instead of Ubuntu's:
<https://docs.docker.com/engine/install/ubuntu/>

On 24.04 this does not arise; `docker.io` there is 29.x.

If you want to run `docker` without `sudo`, add yourself to the group
and start a new login shell — but know that **membership of the `docker`
group is equivalent to root**, since anything that can reach the daemon
can ask it for a privileged container mounting the host filesystem:

```bash
sudo usermod -aG docker "$USER"
```

### 2. Pre-flight (2 minutes, saves an hour)

The base image tag is pinned in the Dockerfile but is worth confirming
against what your host can actually pull, along with the R and
Bioconductor versions it carries:

```bash
docker run --rm bioconductor/bioconductor_docker:RELEASE_3_20 \
  R -q -e 'cat(R.version.string, "| Bioc", as.character(BiocManager::version()), "\n")'
df -h /var/lib/docker   # image is 5-7 GB; build cache wants 2-3x that
```

Then check the pins, which costs a minute and has already paid for
itself twice:

```bash
Rscript deploy/scripts/check_pins.R
```

It resolves the whole dependency closure against the pinned CRAN
snapshot and the Bioconductor mirror and reports any version
requirement that cannot be satisfied. "Every package is present" and
"every package's requirements are satisfiable" are different questions,
and only the second one predicts whether the build works — a snapshot
where all 228 packages existed still died 20 minutes in, on
`BiocParallel` wanting a newer `BH`. Run it after changing
`CRAN_SNAPSHOT`, `BIOC_MIRROR`, the Bioconductor release, or either
package list.

### 3. Build (30-60 minutes)

```bash
deploy/scripts/build_image.sh omicsapp:1.0
```

A failure part way through is not a restart: Docker caches each layer,
so a fix re-runs only from the layer that failed. The Bioconductor and
LaTeX layers are the likely ones. Dropping the LaTeX layer costs the
PDF report and saves 1.5 GB.

### 4. Smoke-test the container on its own — do not skip

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

### 5. Network and storage

```bash
docker network create sp-net
sudo mkdir -p /srv/omicsapp/users            # /srv is the HDD
sudo deploy/scripts/add_user.sh puweilin     # once per user
```

Confirm `/srv` really is the HDD before creating anything — see
[Storage layout](#storage-layout). Creating the directories on the
system disk by mistake is easy to miss and awkward to undo once people
have data in them.

```bash
df -h /srv    # should show the 60 TB volume, not the root filesystem
```

### 6. ShinyProxy

```bash
sudo mkdir -p /etc/shinyproxy
sudo cp deploy/shinyproxy/application.yml.template /etc/shinyproxy/application.yml
sudo chmod 600 /etc/shinyproxy/application.yml
# Owner as well as mode. The .deb's unit runs the service as
# `shinyproxy`, so a root-owned 600 file is one it cannot read: it
# starts, fails on "application.yml (Permission denied)" buried under a
# Spring stack trace, and exits 0 -- which systemd reports as
# "Deactivated successfully".
sudo chown shinyproxy:shinyproxy /etc/shinyproxy/application.yml
```

No passwords go in this file. It names Keycloak's endpoints and carries
one client secret, which you paste in after starting Keycloak — the
order and the reasoning are in [`keycloak/README.md`](keycloak/README.md).

Stand Keycloak up first. ShinyProxy fetches the signing keys at startup
and will not start without them.

**The template ships 8081 and 9091, not the defaults.** This is a shared
machine: 8080 belongs to a colleague's service (`Server: SkinAtlas/1.0`)
and 9090 to something called `mihomo`. ShinyProxy binds two ports — the
one users reach, and a second for Spring's actuator endpoints — and both
conventional choices are taken here.

Check before starting, because neither collision announces itself
usefully:

```bash
sudo ss -tlnp | grep -E ':(8081|9091) ' && echo "TAKEN -- pick another"
```

The actuator collision is the confusing one. The log contradicts itself:
"Undertow started on port 8081" followed by "Web server failed to start.
Port 9090 was already in use." The main port really did come up; the
context failed anyway, so nothing serves.

The wrong-port collision is worse, because it does not fail at all.
Point nginx at 8080 here and it forwards to the neighbouring service,
which answers `401` — a login prompt for software you have never heard
of, and nothing in any log to say ShinyProxy was never involved.

`proxy.port` in `application.yml` and `proxy_pass` in the nginx site are
two files that must agree. After any change to either:

```bash
curl -sI http://127.0.0.1:8081 | head -3
```

If the `Server` header names something else, that port belongs to
another service too.

### 7. nginx and firewall

```bash
sudo cp deploy/nginx/omicsapp.conf /etc/nginx/sites-available/omicsapp
sudo ln -s /etc/nginx/sites-available/omicsapp /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx
sudo ufw allow from 192.168.51.0/24 to any port 80 proto tcp
```

Leave `sites-enabled/default` alone. nginx picks a server block by
matching the request's `Host` against `server_name`, and an exact match
beats `default_server`, so `http://192.168.51.52` reaches this app while
everything else still reaches whatever was there before. Removing the
default site is not needed and takes down any static content a
colleague was serving from it.

Users then reach the app at `http://192.168.51.52`.

### 8. Acceptance

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

The first three below were all latent for the same reason, and it is
worth naming: nobody had ever logged in successfully, so no code path
downstream of "start a container" had ever run. Three separate faults
were sitting in a row, and each one only became visible after the one
before it was fixed. If a deployment has never served a real session,
assume the same and fix them in order rather than concluding the first
one was the problem.

**ShinyProxy before 3.2 cannot start any container on Docker 28+.** The
failure is `Cannot build ImageInfo, some of required attributes are not
set [comment, dockerVersion, author]`, on a request that returned HTTP
200 — the call succeeded and the *parse* failed. Docker no longer
returns those legacy v1 image-config fields, and the `docker-client
7.0.8-OA-3` bundled in 3.1.1 requires them. It affects every image,
pulled or built, so neither rebuilding nor switching off the containerd
image store helps. 3.2.4 bundles `7.0.8-OA-5`, where the three fields
are no longer in the builder's required set. That upgrade also moves
Java 17 to 21.

**`internal-networking` must be false when ShinyProxy runs on the
host.** True makes it address containers by their Docker-network
hostname, which resolves through Docker's embedded DNS — available
inside containers, not to a systemd service. The container starts, the
health check fails with `UnknownHostException: <container id>`, and the
browser says "Failed to start app", which reads like the application
crashed.

**`server.forward-headers-strategy` is required behind TLS.** Without
it every absolute URL is built from the request as nginx forwarded it,
i.e. plain HTTP. nginx redirects those back, so the site works and
merely costs a round trip — but the OIDC `redirect_uri` is built the
same way, and Keycloak rejects a `http://` one against a `https://`
registration.

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

**`customHeader` authentication is available but not used.** It arrived
in 3.2.0, so the pinned 3.2.4 has it, and pairing it with a
forward-auth proxy is a simpler integration than OIDC — no client
secret, no redirect URI, no signing keys to fetch.

It is not used because it makes ShinyProxy trust an unsigned header
absolutely: whoever the header names is who you are. That is safe only
if nothing but the reverse proxy can reach ShinyProxy, and this is a
shared machine — a colleague's service already holds port 8080, so a
local process sending `Remote-User: <someone>` to 127.0.0.1:8081 is a
path that exists here rather than a theoretical one. `bind-address:
127.0.0.1` stops remote hosts, not local ones.

OIDC has no equivalent boundary: the identity is carried in a token
Keycloak signed, so forging it needs Keycloak's key rather than the
ability to open a socket.

**Log in once before adding everybody.** Whatever the doubt is --
version syntax, indentation, the password format -- ten accounts made
wrong fail the same way as one, and one is faster to look at.

## Storage layout

**All user data lives on the HDD, mounted at `/srv`. Docker's own
storage stays on the SSD.**

| What | Path | Disk | Why |
|---|---|---|---|
| Projects and archived uploads | `/srv/omicsapp/users/<user>/` | **HDD** | Grows forever; survives an OS reinstall untouched |
| Docker images and layers | `/var/lib/docker` | **SSD** | Read on every container start; a build is heavy random I/O |
| Backups | `/backup/omicsapp/` | **SSD** | A different physical disk from the data |

Mount by UUID rather than device name — `/dev/sdb` is not stable across
reboots, let alone across a reinstall:

```bash
sudo blkid /dev/sdX1                       # note the UUID
echo 'UUID=<uuid>  /srv  ext4  defaults  0 2' | sudo tee -a /etc/fstab
sudo mount -a && df -h /srv
```

Nothing else in this directory knows which disk it is on. That is
deliberate: a path in a config file that encodes a physical disk has to
be edited to move data, and the two drift.

### Why the HDD, given the SSD has room

Capacity is not the reason — 75 GB a year against 4 TB is fifty years
either way. The reason is **blast radius**. A separate physical disk is
untouched by whatever happens to the system disk: a root partition
filling up, a filesystem repair, a reinstall two years from now. On an
SSD partition the same data would instead depend on somebody
remembering not to tick "format".

Backups only mean anything on a *different* physical disk, so both
disks are in use either way. This only decides which one holds the
original.

The cost is measured and small. A `.omp` file is 5.9 MB and writing one
takes 0.06 s on SSD, most of it serialisation rather than I/O; a
spinning disk adds roughly 40–70 ms. Autosave fires five to seven times
across a full workflow, so the total is a few hundred milliseconds
spread over a session that already spends seconds in `run_diff()`.

### What the numbers mean

Measured on a 20k feature by 60 sample workbook: a `.omp` project file
is **5.9 MB**, the raw upload it came from is **20 MB**, one full
analysis produces about **45 MB**. Three analyses a day is **50–75 GB a
year**, of which roughly 15 GB is archived uploads.

On the 7 GB image, because it changes the arithmetic: layers are
**read-only and shared**. Four concurrent containers do not use 4 × 7 GB
— there is one copy on disk, and each container adds only its writable
layer, which stays tiny because all data goes to the bind mount rather
than into the container.

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

The data is on the HDD, so the backup goes to the **SSD** — a different
physical disk, which is the point. There is ample room: 75 GB a year
against 4 TB.

**Two directories, not one.** `users/` is the work; `keycloak-db/` is
the accounts. Without the second, a restore gives you every file and
nobody who can log in to reach it — and the directory names are UUIDs,
so there is no way to work out whose is whose after the fact.

```bash
sudo mkdir -p /backup/omicsapp
sudo rsync -a --delete /srv/omicsapp/users/       /backup/omicsapp/users/
sudo rsync -a --delete /srv/omicsapp/keycloak-db/ /backup/omicsapp/keycloak-db/
```

```cron
# /etc/cron.d/omicsapp-backup — nightly at 02:30
30 2 * * *  root  rsync -a --delete /srv/omicsapp/users/ /backup/omicsapp/users/
35 2 * * *  root  rsync -a --delete /srv/omicsapp/keycloak-db/ /backup/omicsapp/keycloak-db/
```

Copying Postgres' files while it is running gives you a backup that may
need crash recovery to open. That is usually fine and is much better
than nothing, but if these accounts ever become hard to recreate, take a
real dump instead:

```cron
25 2 * * *  root  docker exec keycloak-db pg_dump -U keycloak keycloak | gzip > /backup/omicsapp/keycloak.sql.gz
```

Two things this does not cover, so decide about them explicitly:

* **An OS reinstall wipes the SSD, and the backup with it.** The data
  itself survives on the HDD, which is why it is there — but you would
  be running without a backup until the SSD is repopulated. Copy
  `/backup` somewhere else before reinstalling.
* **`--delete` mirrors deletions.** A file removed on Monday is gone
  from the backup on Tuesday. That is a mirror, not history; if you want
  to recover a project someone deleted last week, use a snapshotting
  filesystem or `rsync --link-dest` rotations instead.

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

## Reinstalling the host

Not planned. It is documented anyway, because what a reinstall *would*
cost is what justifies the storage layout above, and because the
question comes back every time the host OS is discussed.

Almost nothing here depends on the host distribution. The image carries
its own userland — R, 224 packages, the Bioconductor stack — so the host
needs only Docker, nginx, rsync and `useradd`, which are the same on any
recent Ubuntu.

The container's own operating system is fixed by the `FROM` line, not by
the host: `bioconductor/bioconductor_docker:RELEASE_3_20` is **Ubuntu
24.04, with R 4.4.2 and Bioconductor 3.20**, whatever the host runs. An
image built on one Ubuntu runs on another, and rebuilding on a different
host produces the identical container — so a host change is never by
itself a reason to rebuild.

The host's one real contribution is the Docker engine. On a host older
than 24.04 check its version (step 1): a 24.04 userland on an engine
older than 20.10.10 hits the `clone3` seccomp block, and the symptom is
a container that dies on start rather than an error that says so.

**Export the image before touching the system disk.** Building is the
one step with real unknowns — an hour, and the first attempt usually
turns up a package or version problem. The export is minutes and the
result is portable, so there is no reason to pay the hour twice:

```bash
sudo mkdir -p /srv/backup && sudo chown "$USER" /srv/backup
sudo docker save omicsapp:1.0 | gzip > /srv/backup/omicsapp-1.0.tar.gz
df -h /srv/backup && ls -lh /srv/backup   # on the HDD, not on /

# afterwards
gunzip -c /srv/backup/omicsapp-1.0.tar.gz | sudo docker load
```

Worth doing once the deployment is accepted, reinstall or not: it is the
only copy of the image that is not inside `/var/lib/docker`.

The redirect is deliberate: `>` is performed by *your* shell before
`sudo` runs, so the archive is written as you, not as root. That is why
the directory is chowned first — `sudo docker save > /some/root/path`
fails with a permission error that appears to come from `sudo` and does
not.

What to preserve, in order of how much it hurts to lose:

| | Where | Note |
|---|---|---|
| User data | `/srv/omicsapp/users` | On the HDD; do not format that disk |
| ShinyProxy config | `/etc/shinyproxy/application.yml` | Contains password hashes — copy with mode 600 |
| The image | `docker save` | Rebuildable, but that is an hour |
| nginx site | `/etc/nginx/sites-available/omicsapp` | Also in this repository |
| The code | — | On GitHub; nothing to do |

Docker itself is **not** on that list. It lives on the system disk and
goes with it, along with everything in `/var/lib/docker` — which is
exactly why the image has to be exported to the HDD rather than left
where Docker keeps it. Reinstalling the package is five minutes; the
7 GB it used to hold is the hour.

**The bigger risk is not this application.** Other people's work lives
on that machine — conda environments, half-finished jobs, data that was
never anywhere else. Ask each of them before the disk is touched;
recovering our stack afterwards is an afternoon, recovering theirs may
be impossible.

After the reinstall the system disk is empty, so Docker, nginx and
ShinyProxy all have to go back on — step 1 in full. `/etc/fstab` needs
the HDD entry again by UUID, and the cron backup needs re-adding.

What you skip is the expensive part: `docker load` replaces steps 2 and
3, and the data is already sitting on `/srv`. Everything else is steps
5 to 8, which is about half an hour.

## Hardening, once it works

Leave these until the acceptance checklist passes. First deployments go
wrong for ordinary reasons, and every extra variable is one more
candidate.

**ShinyProxy already has its own account.** The `.deb` creates a
`shinyproxy` system user, adds it to `docker`, and ships a unit that
uses it:

```
User=shinyproxy   Group=shinyproxy   WorkingDirectory=/etc/shinyproxy
uid=997(shinyproxy) gid=1004(shinyproxy) groups=1004(shinyproxy),125(docker)
```

So there is nothing to do here, only something to know — and one thing
to get right, which is why step 6 chowns the config: a root-owned 600
`application.yml` is one the service cannot read.

A dedicated account rather than a person's is the right shape anyway: a
stolen SSH key should not hand over the service, and a compromised
service should not reach someone's files.

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
