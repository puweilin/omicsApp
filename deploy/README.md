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

```bash
# 1. Build the image (30-60 min the first time)
deploy/scripts/build_image.sh omicsapp:1.0

# 2. Network and storage
docker network create sp-net
sudo mkdir -p /srv/omicsapp/users
sudo deploy/scripts/add_user.sh puweilin

# 3. ShinyProxy config
sudo mkdir -p /etc/shinyproxy
sudo cp deploy/shinyproxy/application.yml.template /etc/shinyproxy/application.yml
sudo chmod 600 /etc/shinyproxy/application.yml
# ...fill in the bcrypt password hashes, then start ShinyProxy

# 4. nginx
sudo cp deploy/nginx/omicsapp.conf /etc/nginx/sites-available/omicsapp
sudo ln -s /etc/nginx/sites-available/omicsapp /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx

# 5. Firewall — restrict to the local subnet rather than opening the port
sudo ufw allow from 192.168.51.0/24 to any port 80 proto tcp
```

Users then reach the app at `http://192.168.51.52`.

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
