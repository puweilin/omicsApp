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
| `host.env.template` | The server's address, the one thing every deployment sets differently. Copy to `host.env` (gitignored) and fill in. |
| `scripts/render.sh` | Writes the four files that carry that address from their `.template` neighbours. |
| `rsync.exclude` | What step 0's `rsync --delete` must leave alone on the server: the secrets, `host.env`, and the rendered files. |
| `shinyproxy/application.yml.template` | Rendered to `application.yml`; copy that to the server and fill in the client secret. |
| `nginx/omicsapp.conf.template` | Reverse proxy: TLS, the WebSocket headers Shiny needs, Keycloak at `/auth`, security headers and rate limits. Rendered to `omicsapp.conf`. |
| `nginx/nginx-limits.conf` | systemd drop-in raising nginx's file-descriptor limit. |
| `keycloak/` | The identity provider: compose file and realm definition (both rendered from `.template`), and its own README. |
| `cron/omicsapp-backup` | Runs the nightly backup and the weekly restore drill. |
| `scripts/backup.sh` | Dated, checksummed snapshots of the work, the account database and the configuration; copied off the host; alerts on failure. |
| `scripts/restore_check.sh` | Restore drill: backup age, checksums, the account dump restored into a scratch database, a project opened. |
| `backup.env.template` | Backup settings (remote host, alert hook, retention). Copy to `/etc/omicsapp/backup.env`. |
| `scripts/build_image.sh` | Builds the image from the right context and tags it `<version>-<commit>`, a name that is never reused. |
| `scripts/rollback.sh` | Switches the image ShinyProxy starts and restarts it — for going back, and for going live. |
| `scripts/pin_base_digests.sh` | Resolves the third-party base images to digests and records them in `docker/base-digests.lock`. |
| `docker/base-digests.lock` | The digest each pinned image was resolved to; the contract test holds the Dockerfile and compose file to it. |
| `scripts/egress.sh` | Firewall rules that stop app containers opening any connection. |
| `systemd/omicsapp-egress.service` | Re-applies those rules at boot. |
| `docker/daemon.json` | Docker daemon settings: log rotation and no-new-privileges for every container (host-wide — read [Hardening](#hardening-once-it-works) first). |
| `logrotate/omicsapp` | Rotation for the per-session app logs and the cron jobs' logs. |
| `scripts/add_user.sh` | Creates an account in Keycloak and its storage directory, in that order. |
| `scripts/list_users.sh` | Maps the UUID directory names back to people. |

## The one thing every deployment sets

Your server's address, written once:

```bash
cp deploy/host.env.template deploy/host.env
$EDITOR deploy/host.env             # OMICSAPP_HOST=<server-ip or DNS name>
deploy/scripts/render.sh
```

Four files carry that address and they have to agree — nginx's
`server_name` and the certificate it expects, Keycloak's `KC_HOSTNAME`,
the realm's redirect URIs, ShinyProxy's OIDC endpoints. Each is tracked
as a `.template` with `@OMICSAPP_HOST@` where the address goes;
`render.sh` writes the real file beside it.

So `deploy/host.env` is the only file you edit that is about *your*
machine rather than about the software. It and the four rendered files
are gitignored and excluded from the sync below: they describe this
server, and the next person's server is not this one.

Re-run `render.sh` after a sync that changed a template, and after any
change of address — then re-copy the nginx site, bring Keycloak up again
with the new compose file, and re-import the realm, because the address
is baked into all three.

`host.env` carries one optional second setting, `ADMIN_ALLOW_CIDR`: the
addresses allowed to reach the administrative pages through nginx
(Keycloak's admin console and master realm, ShinyProxy's `/admin`).
Empty, it is the server itself — `127.0.0.1`, `::1` and `OMICSAPP_HOST`.
See [The admin console is not on the LAN](#the-admin-console-is-not-on-the-lan).
A third, `RATE_LIMIT_EXEMPT_CIDR`, lists the addresses nginx's rate
limits do not count, with the same default; see
[Security headers and rate limits](#security-headers-and-rate-limits).

Everything else below is the same on any machine.

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
git rev-parse HEAD > REVISION   # the commit, for the image tag (gitignored)
```

**Getting it onto the server.** There are no git credentials there, so
this is a copy from the machine you work on, not a `git pull`:

```bash
rsync -avz --delete --exclude-from=deploy/rsync.exclude \
  ~/SciProject/CHISSS/omicsApp/ <user>@<server-ip>:~/omicsApp/
```

`REVISION` travels with the copy because `.git` does not: it is how
`build_image.sh` on the server knows which commit it is building, and so
what to call the image (step 3).

Both flags earn their place. `--delete` is what removes a script that
was deleted upstream — otherwise it lingers on the server and someone
runs it a year later. The exclude list protects what exists only on the
server: `deploy/keycloak/.env` with the Keycloak secrets,
`deploy/host.env` with the address, and the four files `render.sh`
writes from them (step 1). rsync does not delete excluded paths, so
naming them there is what keeps them.

If the copy is slow, that is the host's network link, not rsync or the
laptop; `-z` is in the command because the source is text and compresses
several times over.

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

**The server's address** is set once, at the top of this file, before
any of these steps. If `render.sh` has not been run, the files the later
steps copy still say `@OMICSAPP_HOST@` and nothing will resolve.

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

The base image is pinned by **digest**, not only by tag: `ARG
BASE_IMAGE=bioconductor/bioconductor_docker:RELEASE_3_20@sha256:...` in
the Dockerfile, and `postgres:16@sha256:...` in the Keycloak compose
file. A tag is a name its owner can move; a digest names the bytes, so a
rebuild next year starts from the image this one was tested on.
`deploy/scripts/pin_base_digests.sh` asks the registries and rewrites
both, with the record in `deploy/docker/base-digests.lock`;
`test-deploy-contract.R` fails if a reference and the lock disagree, and
warns about any image with no recorded digest. Re-run it when you
*choose* to take an update — a Postgres security release, a rebuilt
base — then rebuild and test as for any change.

> **TODO — Keycloak's digest.** `quay.io/keycloak/keycloak:26.7.3` is
> still pinned by tag only: quay.io could not be reached when the others
> were recorded. Run `deploy/scripts/pin_base_digests.sh` once from a
> machine that can reach quay.io and commit the result; the contract
> test warns until then.

### 3. Build (30-60 minutes)

```bash
deploy/scripts/build_image.sh
# -> Built omicsapp:0.0.0.9000-1a2b3c4d5e6f
```

The tag is `<omicsApp version>-<first 12 characters of the commit>`, and
it is never reused: building a commit that already has an image stops
rather than overwrite it, and uncommitted changes get a `-dirty.<time>`
suffix so an image that matches no commit cannot pass for one that
does. That is what makes [rolling back](#rolling-back) possible —
yesterday's image is still on the host under its own name. Extra
arguments add aliases (`build_image.sh omicsapp:1.0`); `latest` is
refused, because it names whatever was built last.

A failure part way through is not a restart: Docker caches each layer,
so a fix re-runs only from the layer that failed. The Bioconductor and
LaTeX layers are the likely ones. Dropping the LaTeX layer costs the
PDF report and saves 1.5 GB.

### 4. Smoke-test the container on its own — do not skip

This separates "does the app work" from "is ShinyProxy configured
right". Debugging both at once is what makes a deployment take a day.

```bash
TAG=omicsapp:0.0.0.9000-1a2b3c4d5e6f      # what build_image.sh printed
docker run --rm -p 3838:3838 -e OMICSAPP_DATA_DIR=/tmp/data "$TAG"
```

Open `http://<host>:3838` and confirm all seven views render and the
Report view shows an "Analysis code" card. Then check the two things
that fail silently rather than loudly:

```bash
# Gene-set cache is in the image. Missing, the app still works -- it
# just pays ~10s on every user's first enrichment, forever.
docker run --rm "$TAG" \
  R -q -e 'cat(length(list.files(Sys.getenv("OMICSCORE_GENESET_CACHE"))), "cached tables\n")'
# expect 14

# Whether forked workers are available. launch() reads exactly this to
# decide between multicore and multisession (see run_app.R); FALSE costs
# a process spawn per analysis rather than correctness, but it is worth
# knowing which one the container will pick.
#
# NOT `class(future::plan())`: the plan is set by launch(), so a bare
# `R -e` session that never calls it reports `sequential` whatever the
# image contains. That check can only ever fail, which makes it a check
# people learn to ignore.
docker run --rm "$TAG" \
  R -q -e 'cat("supportsMulticore:", future::supportsMulticore(), "\n")'
# expect TRUE
```

### 5. Network and storage

```bash
# Two networks, deliberately. Keycloak and its Postgres share sp-net;
# app containers get their own. An app container has no reason to reach
# the account database, and on one shared network it could -- code
# running user-supplied data would be a hostname away from the
# credentials. Neither container needs Docker DNS for anything else:
# ShinyProxy reaches apps through a published port.
docker network create sp-net
# A fixed bridge name, so the egress rules below can name it and still
# match after the network is recreated.
docker network create -o com.docker.network.bridge.name=br-omicsapp omicsapp-net

# App containers answer ShinyProxy and open nothing themselves -- not the
# internet, not the LAN, not Keycloak or its database, not the host's own
# services. See "Hardening" for why, and how to check.
sudo deploy/scripts/egress.sh apply
sudo cp deploy/systemd/omicsapp-egress.service /etc/systemd/system/
sudo systemctl daemon-reload && sudo systemctl enable --now omicsapp-egress

sudo mkdir -p /srv/omicsapp/users            # /srv is the HDD
sudo deploy/scripts/add_user.sh you@example.com   # once per person
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
sudo cp deploy/shinyproxy/application.yml /etc/shinyproxy/application.yml   # rendered in step 1
sudo chmod 600 /etc/shinyproxy/application.yml
# Owner as well as mode. The .deb's unit runs the service as
# `shinyproxy`, so a root-owned 600 file is one it cannot read: it
# starts, fails on "application.yml (Permission denied)" buried under a
# Spring stack trace, and exits 0 -- which systemd reports as
# "Deactivated successfully".
sudo chown shinyproxy:shinyproxy /etc/shinyproxy/application.yml

# Where ShinyProxy writes each app session's output (container-log-path),
# and its rotation.
sudo install -d -o shinyproxy -g shinyproxy -m 750 /var/log/shinyproxy/containers
sudo cp deploy/logrotate/omicsapp /etc/logrotate.d/omicsapp

# The image users get: the tag build_image.sh printed (step 3).
sudo sed -i -E "s|^([[:space:]]*container-image:[[:space:]]*).*|\1$TAG|" /etc/shinyproxy/application.yml
```

Later image changes — upgrades and rollbacks alike — go through
`deploy/scripts/rollback.sh`, which also restarts ShinyProxy; this once
it is not running yet.

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

# Raise the file-descriptor limit before the first reload, not after the
# first failure -- see the file's own comments for what that looks like.
sudo mkdir -p /etc/systemd/system/nginx.service.d
sudo cp deploy/nginx/nginx-limits.conf /etc/systemd/system/nginx.service.d/limits.conf
sudo systemctl daemon-reload

sudo nginx -t && sudo systemctl restart nginx
sudo ufw allow from <lan-cidr> to any port 443 proto tcp
sudo ufw allow from <lan-cidr> to any port 80 proto tcp
```

`restart`, not `reload`, for that first one: a reload keeps the running
master, and the master is what holds the old limit.

Then check the admin pages are closed to the LAN — from your
workstation, not the server:

```bash
curl -sk -o /dev/null -w '%{http_code}\n' https://<server-ip>/auth/admin/          # 403
curl -sk -o /dev/null -w '%{http_code}\n' https://<server-ip>/auth/realms/master/  # 403
curl -sk -o /dev/null -w '%{http_code}\n' https://<server-ip>/auth/realms/omicsapp/account  # 200 or 302
```

And that the security headers arrive, once each:

```bash
curl -skI https://<server-ip>/ | grep -iE 'strict-transport|x-frame|x-content-type|referrer|permissions|content-security'
# Strict-Transport-Security: max-age=31536000
# X-Content-Type-Options: nosniff
# X-Frame-Options: SAMEORIGIN
# Referrer-Policy: same-origin
# Permissions-Policy: camera=(), ...
# Content-Security-Policy: frame-ancestors 'self'; object-src 'none'; base-uri 'self'
# Content-Security-Policy-Report-Only: default-src 'self'; ...
```

`<lan-cidr>` is the network's real prefix, read from `ip -br addr` on
the host rather than assumed. A rule written as a `/24` on a `/21`
admits a quarter of the building and refuses the rest, and it stays
latent for exactly as long as ufw is not enforcing.

Confirm it took:

```bash
grep 'Max open files' /proc/$(cat /var/run/nginx.pid)/limits
```

Leave `sites-enabled/default` alone. nginx picks a server block by
matching the request's `Host` against `server_name`, and an exact match
beats `default_server`, so `http://<server-ip>` reaches this app while
everything else still reaches whatever was there before. Removing the
default site is not needed and takes down any static content a
colleague was serving from it.

Users then reach the app at `https://<server-ip>`; plain `http://` redirects there.

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

**A slow copy to the host is the host's link, not rsync.** When the
step 0 sync crawls, measure the direction before blaming the tool:
inbound to a host whose only link is Wi-Fi can be an order of magnitude
slower than outbound from it, whoever sends. The measurements for this
deployment, and what changing the link would do to the address the four
rendered files carry, are in the lab's internal notes rather than here.

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
| `OMICSAPP_SIGNING_KEY` | unset | Key the store signs project files with. Unset, each user's store generates one and keeps it as `.omicsapp-signing-key`. See [Project files are signed](#project-files-are-signed). |
| `OMICSAPP_MAX_PROCS` | `2048` | Soft cap on processes for the app's uid (`ulimit -u`), standing in for `--pids-limit`. Per uid, shared by every user's container. |

## Operations

**Backup — do this on day one.** Losing a disk is far likelier than any
attack, and `.omp` files are the irreplaceable part: a raw upload can
usually be exported from the instrument again, the analysis parameters
encoded in a project cannot.

`deploy/scripts/backup.sh` does the whole job; the header of the script
says what it keeps and why. In short, every night it writes a **dated
snapshot** of:

* `/srv/omicsapp/users/` — the work;
* the Keycloak account database, as a `pg_dump` that is checked before it
  replaces anything (and the database files as a fallback). Without it a
  restore gives you every file and nobody who can log in to reach it —
  the directory names are UUIDs, so there is no way to work out whose is
  whose after the fact;
* the rendered configuration (`application.yml`, nginx, Keycloak `.env`
  and realm, `host.env`), the TLS certificate and key, the cron files and
  the refreshed gene-set cache — what it takes to stand the service up
  again, plus `/etc/docker/daemon.json`, the logrotate file and the
  egress unit;
* the logs (`LOG_PATHS`): ShinyProxy's own and every app session's
  (`/var/log/shinyproxy`), nginx's, and these jobs' — so the record of
  what happened survives the disk it was written on;
* `MANIFEST.sha256` (every file's checksum) and `IMAGE` (the digest and
  tag of the image ShinyProxy is configured to run, read from
  `application.yml`).

Each user's directory also holds the store's signing key
(`.omicsapp-signing-key`) when `OMICSAPP_SIGNING_KEY` is not set, so the
snapshot of `users/` carries what is needed to open its projects.

Unchanged files are hard links to the previous snapshot, so a night costs
only what changed. Snapshots are kept for 14 days, then one per week for
8 weeks. With `BACKUP_REMOTE` set they are also copied, with their
history, to **another machine** — without it the backup shares the
host's fate (fire, theft, an OS reinstall). Any failure is sent to
`ALERT_WEBHOOK` and/or `ALERT_EMAIL` and logged; the script exits
non-zero and the earlier snapshots are untouched.

```bash
sudo mkdir -p /etc/omicsapp /backup/omicsapp
sudo cp deploy/backup.env.template /etc/omicsapp/backup.env
sudo chmod 600 /etc/omicsapp/backup.env
sudo $EDITOR /etc/omicsapp/backup.env       # BACKUP_REMOTE, ALERT_WEBHOOK, ...
sudo deploy/scripts/backup.sh               # once by hand: it should end with "done"
sudo cp deploy/cron/omicsapp-backup /etc/cron.d/omicsapp-backup
sudo chmod 644 /etc/cron.d/omicsapp-backup
```

The cron file expects the deployed checkout at `/opt/omicsApp` (set
`REPO_DIR` and the paths in the cron file if it is elsewhere). For
`BACKUP_REMOTE`, give root an ssh key the other machine accepts — cron
cannot answer a password prompt.

`chmod 644` is not decoration: cron silently ignores a file in
`/etc/cron.d` that is group- or world-writable, so a stricter-looking
mode gets you no backups and no error.

**The restore drill.** A backup that has never been restored is a hope.
Every Sunday the cron file runs `deploy/scripts/restore_check.sh`, which
alerts when the newest backup is more than 36 hours old, verifies the
checksums, restores the Keycloak dump into a throwaway Postgres container
and counts the accounts, opens a random project with the image's own
`load_project()`, and checks that the remote copy has the latest snapshot.
Run it by hand after setting up, and after any change to the server:

```bash
sudo deploy/scripts/restore_check.sh        # ends with "all checks passed"
```

**Restoring.** Stop ShinyProxy and Keycloak, then from the snapshot you
want (`/backup/omicsapp/latest` or a dated one under `snapshots/`):

```bash
S=/backup/omicsapp/latest
(cd "$S" && sha256sum --quiet -c MANIFEST.sha256)          # nothing printed = intact
sudo rsync -a "$S/users/" /srv/omicsapp/users/
docker compose -f deploy/keycloak/docker-compose.yml up -d keycloak-db
gunzip -c "$S/keycloak.sql.gz" | docker exec -i keycloak-db psql -U keycloak keycloak
# configuration, if the host itself was lost:
ls "$S/config"     # deploy/... files and host/etc/... (certificate, cron)
```

Hard-linked snapshots share unchanged files, so never edit a file inside
a snapshot: copy it out first.

**Keeping KEGG current.** The image bakes the MSigDB tables in, and the
`kegg` table among them is the 2011 `KEGG_LEGACY` snapshot (186 human
pathways; KEGG today has ~370). `omicsCore::refresh_geneset_cache()`
replaces it with a live fetch from the KEGG REST API and re-snapshots
the other databases from msigdbr. To keep that current without
rebuilding images, move the cache onto a volume and refresh it on a
schedule:

```bash
# one-time: seed the volume, then mount it in application.yml with
#   container-volumes: [ ..., "/srv/omicsapp/genesets:/opt/genesets:ro" ]
mkdir -p /srv/omicsapp/genesets
docker run --rm -v /srv/omicsapp/genesets:/opt/genesets omicsapp:1.0 \
  R -q -e 'for (org in c("Hs","Mm")) omicsCore::refresh_geneset_cache(organism = org, force = TRUE)'
```

```bash
# /etc/cron.d/omicsapp-genesets — monthly, 03:00 on the 1st
0 3 1 * * root docker run --rm -v /srv/omicsapp/genesets:/opt/genesets omicsapp:1.0 R -q -e 'for (org in c("Hs","Mm")) omicsCore::refresh_geneset_cache(organism = org, force = TRUE)' >> /var/log/omicsapp-genesets.log 2>&1
```

(Use the tag ShinyProxy runs — `grep container-image
/etc/shinyproxy/application.yml` — in place of `omicsapp:1.0`.) These
runs are on Docker's default network, not `omicsapp-net`, which is why
they can reach KEGG while the app containers cannot; the app's
`OMICSCORE_GENESET_TTL_DAYS: "0"` in `application.yml` stops it trying.

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
git rev-parse HEAD > REVISION && rsync ...      # on your machine, as in step 0
deploy/scripts/build_image.sh                   # on the server
sudo deploy/scripts/rollback.sh <the tag it printed>
```

### Rolling back

Every build keeps its own tag, `<version>-<commit>`, and `rollback.sh`
puts any of them live:

```bash
sudo deploy/scripts/rollback.sh --list     # images on this host, the one running, recent switches
sudo deploy/scripts/rollback.sh omicsapp:0.0.0.9000-1a2b3c4d5e6f
```

It refuses an image that is not on the host, copies `application.yml`
to `application.yml.bak-<time>`, rewrites its `container-image:` line,
appends `<time> <old> -> <new>` to `/etc/shinyproxy/image-history`,
restarts ShinyProxy and waits for it to answer. Restarting ends every
running session — the app autosaves, so users reload and restore — so
it asks first; `--yes` skips the question. Going live with a new build
is the same command, which is why the old tag is always one
`--list` away.

What keeps the old images there: nothing deletes them but you. Do not
run `docker image prune -a` on this host — between sessions no
container uses the app image, so it counts as unused and goes.
`docker save` the ones you want to keep beyond the disk (see
[Reinstalling the host](#reinstalling-the-host)); `rollback.sh` tells
you to `docker load` one that is missing.

Projects need no rollback of their own. A file written by a newer
release opens in an older one: the signature trailer is invisible to a
release that predates it, and `schema_version` stays at `1.x` for as
long as older readers can cope (each `.omp` also records an integer
`format_version`, and `load_project()` upgrades older files through
registered migrations). One wrinkle, rolling *forward* again past the
release that introduced signing: projects the older release saved in
the meantime are unsigned, and the store now refuses them. Let it adopt
them by re-running its one-time check (next section):

```bash
sudo rm /srv/omicsapp/users/*/.omicsapp-signatures
```

### Project files are signed

The app opens only project files it wrote. A `.omp` is a serialised R
object, and reading one builds whatever it describes, so a file put in
someone's directory by anything other than the app — a copied-in file, a
restore from a doubtful source, another process on the host — is a way
to run code as that user. Every project and autosave the store writes
ends in an HMAC-SHA256 signature, which is checked *before* the file is
read; an unsigned or altered file is refused with "was not saved by this
app" or "has been changed since this app saved it".

**The key.** Set `OMICSAPP_SIGNING_KEY` under `container-env` in
`application.yml` (the commented line is there), on the server only,
like the client secret:

```bash
openssl rand -hex 32
```

Unset, each user's store generates its own key on first use and keeps
it in `/srv/omicsapp/users/<sub>/.omicsapp-signing-key` (mode 600). That
works, but the key then sits beside the files it protects, so anything
that can write there can read it too; the setting is the stronger
arrangement. Switching to it later is safe: on its next use a store
re-signs its files under the new key and deletes the old key file.
Removing the setting again is not — the stores would generate fresh keys
and refuse their own files until it is put back.

**Existing stores.** The first time a store is used by a release with
signing, it signs every existing unsigned project and snapshot that
passes a structure check — data only, no functions, environments,
external pointers or unexpected classes — and writes
`.omicsapp-signatures` to record that the check ran. The project bytes
are not rewritten, only signed. A file that fails the check is left
alone, refused when opened, and named in the container log
(`/var/log/shinyproxy/containers/`). After that one run, an unsigned
file is refused: the store did not write it. Delete
`.omicsapp-signatures` to run the check again, e.g. after a rollback.

There is no other way into the app for a project file: it opens projects
only from the user's own store, never from an upload.

### Where the logs are

All on the host, all rotated, all in the backup:

| What | Where | Rotated by |
|---|---|---|
| ShinyProxy | `/var/log/shinyproxy/shinyproxy.log` | Spring Boot itself (`logging.logback.rollingpolicy` in `application.yml`): 10 MB files, 30 days, 1 GB cap |
| Each app session's R output | `/var/log/shinyproxy/containers/` (`proxy.container-log-path`) | `deploy/logrotate/omicsapp` |
| nginx | `/var/log/nginx/omicsapp.{access,error}.log` | the distribution's logrotate |
| Keycloak and Postgres | `docker logs keycloak` / `keycloak-db` | the compose file's `logging:` (local driver, 5 × 10 MB) |
| Backup and restore drill, gene-set refresh | `/var/log/omicsapp-backup.log`, `/var/log/omicsapp-genesets.log` | `deploy/logrotate/omicsapp` |

The per-session files matter more than they look: ShinyProxy removes a
container when its session ends, and `docker logs` goes with it, so
without `container-log-path` the R error behind "it crashed yesterday"
no longer exists anywhere.

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
| ShinyProxy config | `/etc/shinyproxy/application.yml` | Contains the client secret and, if set, `OMICSAPP_SIGNING_KEY` — copy with mode 600. Without that key the projects will not open |
| The image | `docker save` | Rebuildable, but that is an hour |
| nginx site | `/etc/nginx/sites-available/omicsapp` | Rendered from `nginx/omicsapp.conf.template` and `host.env` |
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

### The admin console is not on the LAN

nginx serves Keycloak's admin console and admin API (`/auth/admin/`),
the master realm its administrator signs in through
(`/auth/realms/master/`), Keycloak's health and metrics paths,
ShinyProxy's admin page (`/admin`) and Spring's actuator paths only to
the addresses in `ADMIN_ALLOW_CIDR` (`deploy/host.env`); everyone else
gets 403. Users never need any of them — their own account page,
`/auth/realms/omicsapp/account`, stays open. The check uses nginx's
normalised path, the same one it forwards, so `/auth//admin` or
`/auth/%61dmin` cannot get past it. Keycloak (8180), ShinyProxy (8081)
and its actuator (9091) listen on 127.0.0.1 only, so nginx is the only
way to them from the network.

The default list is the server itself. To use the console:

* **through the server** — `ssh -D 1080 <you>@<server>`, set the
  browser's SOCKS proxy to `localhost:1080`, and open
  `https://<server>/auth/admin/`. The request then comes from the
  server's own address, which is on the list. Nothing to change;
* **from a fixed workstation** — add its address to `host.env`
  (`ADMIN_ALLOW_CIDR="127.0.0.1 ::1 <server-ip> <workstation-ip>"`),
  re-run `render.sh`, re-copy the nginx site and `sudo nginx -t && sudo
  systemctl reload nginx`.

`add_user.sh` and `list_users.sh` talk to Keycloak on 127.0.0.1:8180
directly and are unaffected.

### Security headers and rate limits

nginx adds these to every HTTPS response, error pages included
(`nginx/omicsapp.conf.template` explains each in place):

| Header | Value | Why |
|---|---|---|
| `Strict-Transport-Security` | `max-age=31536000` | After one visit the browser goes straight to `https://` for a year. No `includeSubDomains` (the server does not own its parent domain's other names) and no `preload` (months to undo). Browsers ignore it on a bare IP address and over a certificate they were told to accept, so on a LAN with a self-signed certificate it does nothing until a real certificate is in place. |
| `X-Content-Type-Options` | `nosniff` | A file is only ever run as what the server said it is. |
| `X-Frame-Options` | `SAMEORIGIN` | ShinyProxy shows the app in an iframe of its own page on the same origin; no other site may frame it. |
| `Referrer-Policy` | `same-origin` | Links out of the app (a pathway database, say) do not carry the app's internal URLs with them. |
| `Permissions-Policy` | `camera=(), microphone=(), geolocation=(), payment=(), usb=()` | Nothing here uses them. |
| `Content-Security-Policy` | `frame-ancestors 'self'; object-src 'none'; base-uri 'self'` | Enforced: the part known to be safe for Shiny and ShinyProxy. |
| `Content-Security-Policy-Report-Only` | `default-src 'self'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; ...` | **Not enforced** — see below. |

ShinyProxy and Keycloak send some of these themselves, with other
values (both send HSTS with `includeSubDomains`; Spring's default
framing rule is `DENY`). nginx hides theirs and sends its own, so the
browser sees each header once. Under `/auth/` it leaves Keycloak's
`X-Frame-Options`, `Content-Security-Policy` and `Referrer-Policy` alone:
Keycloak tunes those per realm (*Realm settings → Security defenses*)
to what its login pages and admin console need.

**The full Content-Security-Policy is report-only.** Shiny needs inline
scripts and styles, the table and plot widgets evaluate JavaScript from
strings, and ShinyProxy's own pages carry inline script, so a policy
that allows all of that is still not proven to allow everything — and
an enforced policy that misses one thing gives a page that loads and
then does nothing. In report-only mode the browser blocks nothing and
writes what it *would* have blocked to the developer console
(F12 → Console, "Content-Security-Policy-Report-Only"). To enforce it:
go through a full session — sign in, import, run each analysis, open
and download a report, sign out — with the console open; if nothing is
reported, rename `Content-Security-Policy-Report-Only` to
`Content-Security-Policy` in the template (merging it with the
enforced one), re-render, re-copy the site and reload nginx.

Rate limits, per client address, answered with `429` when spent:

| Budget | What | Rate | Burst |
|---|---|---|---|
| `omicsapp_login` | Keycloak's password form (`/auth/realms/*/login-actions/`) and token endpoint | 10 a minute | 10 |
| `omicsapp_general` | every other request | 20 a second | 200, no delay |

A person signing in uses two or three login requests; opening the app
fetches a few dozen files at once, which the burst absorbs. Keycloak's
own brute-force lockout (`bruteForceProtected` in the realm) still
applies on top — that one is per account, this one per address.

Two things are never counted:

* **Shiny's WebSocket.** It is one request that then carries the whole
  session, so its traffic is not counted anyway; the upgrade itself is
  exempt too, so a tab reconnecting in a loop after a network blip does
  not spend its user's budget and grey out the app with a 429 nobody
  sees.
* **The server itself.** ShinyProxy fetches every user's tokens from
  Keycloak through nginx (its `token-url` is the public address), from
  the server's own address. Counted, that one address would lock
  everyone out of signing in together. The list is
  `RATE_LIMIT_EXEMPT_CIDR` in `host.env`; empty, it is `127.0.0.1`,
  `::1` and `OMICSAPP_HOST`. If `OMICSAPP_HOST` is a name that resolves
  differently on the server than for users, put the server's real
  address there.

Users behind one NAT share one budget. At this deployment's size that
is still far from the limit; if 429s ever appear in
`/var/log/nginx/omicsapp.error.log` ("limiting requests") during normal
use, raise the burst in the template rather than exempting the NAT,
which would exempt everyone behind it.

### What confines an app container

| | How | Where |
|---|---|---|
| Memory, CPU | `container-memory-limit: 6g`, `container-cpu-limit: 2` | `application.yml` |
| No privileged mode | `container-privileged: false` | `application.yml` |
| Ports on loopback only | `proxy.docker.target-bind-ip: 127.0.0.1` | `application.yml` |
| No outbound connections | `omicsapp-net` + `scripts/egress.sh` (below) | step 5 |
| No capabilities, no privilege gain | runs as an unprivileged user; the image strips every setuid/setgid bit and file capability, and the build fails if one survives | `docker/Dockerfile` |
| Process cap | `ulimit -u ${OMICSAPP_MAX_PROCS:-2048}` before R starts | `docker/Dockerfile` (`CMD`) |
| no-new-privileges, log rotation | daemon-wide defaults | `docker/daemon.json` (below) |
| Read-only root filesystem | not available through ShinyProxy; CI runs the image with `--read-only` and tmpfs `/tmp`, `/home/omics`, `/data` so it is known to work | `.github/workflows/production-image.yaml` |

The gaps are ShinyProxy's, not choices: its Docker backend (3.2.x,
`ContainerSpec`/`DockerEngineBackend` in containerproxy) passes memory
and CPU limits and requests, networks, DNS, volumes, environment,
labels, `privileged`, `docker-ipc`, `docker-runtime`, `docker-user`,
`docker-group-add` and device requests — and nothing for `--cap-drop`,
`--security-opt`, `--read-only`, `--pids-limit` or `--tmpfs`. Keys for
those in `application.yml` would look like hardening and do nothing, so
the table above puts each protection where it can actually take effect.
The process cap is per uid rather than per container (every user's
container runs as uid 1001), which is why it is generous.

The Keycloak containers, which compose starts, get the real thing:
`no-new-privileges`, `cap_drop: [ALL]` (Postgres keeps the five its
entrypoint needs to set up its directory and drop to `postgres`),
`pids_limit`, `mem_limit`, `cpus`, rotated logs, and a read-only root
for Postgres with tmpfs for its socket and `/tmp`. Keycloak's own root
stays writable: `start` re-augments the server when its build options
differ from the image's. Apply with `docker compose up -d` in
`deploy/keycloak/`.

**Docker daemon defaults.** `deploy/docker/daemon.json` sets the `local`
log driver with rotation (5 × 10 MB per container) and
`no-new-privileges` for every new container. It is host-wide: on this
shared machine it applies to colleagues' containers too, and
no-new-privileges breaks any of theirs that relies on `sudo` or another
setuid program inside. Ask first; drop that line if anyone needs it.
Merge it into an existing `/etc/docker/daemon.json` rather than
overwriting, then:

```bash
sudo systemctl restart docker     # stops every container on the host
```

Containers created before the restart keep their old settings until
recreated (`docker compose up -d --force-recreate` for Keycloak).

### App containers cannot reach the network

They need nothing outbound: gene sets come from the image or the
pre-warmed volume, and ShinyProxy is the one that connects *in*. So
`deploy/scripts/egress.sh` drops every new connection a container on
`omicsapp-net` opens — in Docker's `DOCKER-USER` chain (the internet, the
LAN, other Docker networks, Keycloak and its database included) and in
`INPUT` (the host's own services: sshd, nginx, a colleague's :8080).
Replies to ShinyProxy's connections are not new connections, so the app
is unaffected. `omicsapp-egress.service` puts the rules back after a
reboot.

Not `docker network create --internal`, which would be simpler: an
internal network has no published ports, and a published port is how a
ShinyProxy running on the host reaches its containers.

Check it, with a session running:

```bash
sudo deploy/scripts/egress.sh check
docker run --rm --network omicsapp-net "$TAG" \
  R -q -e 'tryCatch({readLines(url("https://cran.r-project.org"), 1); cat("OPEN\n")}, error = function(e) cat("blocked\n"))'
# expect: blocked
```

With Docker's opt-in nftables firewall backend there is no
`DOCKER-USER` chain; the script says so, and the equivalent rule goes
into Docker's nftables setup instead.

## What CI checks

Everything in this directory is checked on every push, to any branch,
and on every pull request — not only on `main`. The workflows are in
`.github/workflows/`:

| Workflow | What | Blocking |
|---|---|---|
| `R-CMD-check.yaml` | `R CMD check --as-cran` on Ubuntu with R release and the previous R minor (both packages) and on macOS with R release (omicsCore); then the tests that read the repository rather than the installed package, run from the source tree | yes |
| `production-image.yaml` | Builds `docker/Dockerfile`, smoke-tests it confined, runs both suites inside it (R 4.4.2, Bioconductor 3.20, the pinned CRAN snapshot) | yes |
| `lint.yaml` | `shellcheck` on `scripts/*.sh`; `hadolint` on the Dockerfile (ignored rules and their reasons in `.hadolint.yaml` at the repository root); `actionlint` on the workflows; every action pinned to a commit SHA; `lintr` on both packages | all but `lintr` |
| `coverage.yaml` | omicsCore's test coverage: a per-file table on the run's summary page and `coverage.xml` (Cobertura) as an artifact. No secrets, no threshold | no |
| `nightly.yaml` | 02:37 UTC and on demand: the performance budget (`OMICSCORE_PERF_TESTS=1`) and the heavy-engine fuzz sweep (`OMICSCORE_FUZZ_TESTS=1`) | — |

**The deploy contract runs under check.** `test-deploy-contract.R`
(in omicsApp) reads this directory, and the hygiene tests read the
package sources; R CMD check runs tests from a copy in
`<pkg>.Rcheck/`, where neither is beside them, so they used to skip on
every CI run. CI sets `OMICSAPP_REPO_ROOT` to the checkout and they run;
set to a path without the repository, it is an error rather than a
skip. Unset — an install from a tarball, CRAN-style — they still skip
cleanly. `.github/scripts/run-source-tests.R` runs the remaining
source-reading tests from the source tree and fails if any of them
skipped for want of it. To run the contract locally:

```bash
Rscript -e 'devtools::load_all("packages/omicsApp"); testthat::test_file("packages/omicsApp/tests/testthat/test-deploy-contract.R")'
```

**lintr reports but does not fail yet.** It had never been run over this
code, so its first runs list a backlog, not regressions; the findings
appear as annotations on a pull request's changed lines. `.lintr` at the
repository root allows 120-character lines and turns off the linters
that disagree with the existing style or need the package loaded
(object names, object usage, indentation, commented code). Once the
backlog is cleared, set `LINTR_ERROR_ON_LINT: "true"` in `lint.yaml`.

**Actions are pinned to commit SHAs**, with the release in a comment
(`actions/checkout@<sha> # v4.4.0`): a tag can be moved to other code, a
SHA cannot. Dependabot (`.github/dependabot.yml`) proposes updates
monthly, SHA and comment together. `.github/scripts/pin-actions.sh`
checks that each SHA is what its tag says; `--update` re-pins from the
tags. The linters themselves are fixed versions checked against
recorded SHA-256 sums (`.github/scripts/install-lint-tools.sh`).
