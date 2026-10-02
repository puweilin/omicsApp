#!/usr/bin/env bash
#
# Nightly backup for omicsApp: versioned, verified, copied off the host,
# and loud when it fails.
#
#   sudo cp deploy/backup.env.template /etc/omicsapp/backup.env   # then edit
#   sudo deploy/scripts/backup.sh                                  # run once by hand
#   sudo cp deploy/cron/omicsapp-backup /etc/cron.d/                # then nightly
#
# What the previous cron lines got wrong, and what this does instead:
#
#   * `rsync --delete` into one directory was a mirror, not a backup: a
#     project deleted (or corrupted) on Monday was gone from the backup
#     on Tuesday. Each run is now a dated snapshot; unchanged files are
#     hard links to the previous one, so a snapshot costs only what
#     changed. Daily snapshots are kept for KEEP_DAILY days, and one per
#     week for KEEP_WEEKLY weeks.
#   * `pg_dump | gzip > keycloak-N.sql.gz` without pipefail: when the
#     database container was down, pg_dump failed, gzip succeeded, and a
#     good dump was replaced by an empty one with exit status 0. The
#     dump now goes to a temporary file, is checked, and only then moved
#     into the snapshot.
#   * Backups lived on the same host as the data. With BACKUP_REMOTE set
#     the snapshots are copied to another machine (rsync over ssh), with
#     their history (-H keeps the hard links).
#   * Nobody heard about a failure. Any failure now calls the alert hook
#     (webhook and/or mail) and the script exits non-zero; success writes
#     $BACKUP_ROOT/last_success, which restore_check.sh watches.
#   * Only the user files and the account database were kept. The
#     rendered configuration (application.yml, nginx, Keycloak .env and
#     realm, host.env), the TLS certificate and key, the cron files and
#     the refreshed gene-set cache are part of being able to restore, and
#     are copied too.
#
# Each snapshot carries MANIFEST.sha256 (every file's checksum) and
# IMAGE (the digest of the running app image), so a restore can be
# verified and matched to the code that wrote it.

set -Eeuo pipefail
umask 077

CONF="${OMICSAPP_BACKUP_CONF:-/etc/omicsapp/backup.env}"
# shellcheck disable=SC1090
[ -r "$CONF" ] && . "$CONF"

USERS_DIR="${USERS_DIR:-/srv/omicsapp/users}"
KEYCLOAK_DB_DIR="${KEYCLOAK_DB_DIR:-/srv/omicsapp/keycloak-db}"
GENESETS_DIR="${GENESETS_DIR:-/srv/omicsapp/genesets}"
REPO_DIR="${REPO_DIR:-/opt/omicsApp}"
BACKUP_ROOT="${BACKUP_ROOT:-/backup/omicsapp}"
KEEP_DAILY="${KEEP_DAILY:-14}"
KEEP_WEEKLY="${KEEP_WEEKLY:-8}"
BACKUP_REMOTE="${BACKUP_REMOTE:-}"
KEYCLOAK_CONTAINER="${KEYCLOAK_CONTAINER:-keycloak-db}"
APP_IMAGE="${APP_IMAGE:-omicsapp:1.0}"
ALERT_WEBHOOK="${ALERT_WEBHOOK:-}"
ALERT_EMAIL="${ALERT_EMAIL:-}"
DOCKER="${DOCKER:-docker}"
# Extra files worth keeping, space-separated (absolute paths).
EXTRA_PATHS="${EXTRA_PATHS:-/etc/ssl/certs/omicsapp.crt /etc/ssl/private/omicsapp.key /etc/cron.d/omicsapp-backup /etc/cron.d/omicsapp-genesets}"

SNAP_DIR="$BACKUP_ROOT/snapshots"
STAMP="$(date +%Y-%m-%d_%H%M%S)"
WORK="$SNAP_DIR/$STAMP.partial"
FINAL="$SNAP_DIR/$STAMP"
LATEST="$BACKUP_ROOT/latest"

log() { printf '%s backup: %s\n' "$(date '+%F %T')" "$*" >&2; }

alert() {
  local msg="omicsApp backup on $(hostname): $*"
  log "ALERT: $*"
  command -v logger >/dev/null 2>&1 && logger -t omicsapp-backup -p user.err "$msg" || true
  if [ -n "$ALERT_WEBHOOK" ]; then
    # JSON-escape the two characters that can break the payload.
    local body
    body="$(printf '%s' "$msg" | sed 's/\\/\\\\/g; s/"/\\"/g')"
    curl -fsS -m 20 -H 'Content-Type: application/json' \
      -d "{\"text\":\"$body\"}" "$ALERT_WEBHOOK" >/dev/null 2>&1 || log "webhook failed"
  fi
  if [ -n "$ALERT_EMAIL" ] && command -v mail >/dev/null 2>&1; then
    printf '%s\n' "$msg" | mail -s "omicsApp backup FAILED" "$ALERT_EMAIL" || true
  fi
  if [ -n "${ALERT_LOG:-}" ]; then printf '%s\n' "$msg" >> "$ALERT_LOG"; fi
}

on_error() {
  local status=$? line=$1
  alert "failed at line $line (exit $status); the previous snapshots are untouched"
  rm -rf "$WORK"
  exit "$status"
}
trap 'on_error $LINENO' ERR

mkdir -p "$SNAP_DIR"
# One run at a time: a slow night must not overlap the next one.
exec 9>"$BACKUP_ROOT/.lock"
if ! flock -n 9; then
  alert "another backup is still running; skipped"
  exit 1
fi

link_dest() {
  # rsync --link-dest wants the previous copy of the same subtree.
  if [ -d "$LATEST/$1" ]; then printf -- '--link-dest=%s' "$(readlink -f "$LATEST/$1")"; fi
}

mkdir -p "$WORK"

# ---- 1. the account database, as a dump that restores ------------------
if [ -n "$KEYCLOAK_CONTAINER" ]; then
  tmp="$WORK/.keycloak.sql.gz.tmp"
  "$DOCKER" exec "$KEYCLOAK_CONTAINER" pg_dump -U keycloak keycloak | gzip -c > "$tmp"
  gzip -t "$tmp"
  # An empty database still dumps a few KB of schema; less is a failure
  # pg_dump did not report.
  if [ "$(gzip -cd "$tmp" | head -c 4096 | wc -c)" -lt 512 ]; then
    alert "the Keycloak dump is empty"
    false
  fi
  mv "$tmp" "$WORK/keycloak.sql.gz"
else
  log "KEYCLOAK_CONTAINER is empty: no database dump this run"
fi

# ---- 2. the work, and the database files as a fallback ------------------
for pair in "users:$USERS_DIR" "keycloak-db:$KEYCLOAK_DB_DIR" "genesets:$GENESETS_DIR"; do
  name="${pair%%:*}"
  src="${pair#*:}"
  if [ -d "$src" ]; then
    # shellcheck disable=SC2046
    rsync -a --delete $(link_dest "$name") "$src/" "$WORK/$name/"
  elif [ "$name" = "users" ]; then
    alert "the user directory $src does not exist"
    false
  fi
done

# ---- 3. what it takes to stand the service up again ---------------------
mkdir -p "$WORK/config"
for rel in deploy/host.env deploy/keycloak/.env deploy/keycloak/docker-compose.yml \
           deploy/keycloak/omicsapp-realm.json deploy/shinyproxy/application.yml \
           deploy/nginx/omicsapp.conf; do
  if [ -f "$REPO_DIR/$rel" ]; then
    mkdir -p "$WORK/config/$(dirname "$rel")"
    cp -p "$REPO_DIR/$rel" "$WORK/config/$rel"
  fi
done
for f in $EXTRA_PATHS; do
  if [ -f "$f" ]; then
    mkdir -p "$WORK/config/host$(dirname "$f")"
    cp -p "$f" "$WORK/config/host$f"
  fi
done
if [ -d "$REPO_DIR/.git" ]; then
  git -C "$REPO_DIR" rev-parse HEAD > "$WORK/config/REPO_COMMIT" 2>/dev/null || true
fi
"$DOCKER" image inspect --format '{{.Id}} {{index .RepoTags 0}}' "$APP_IMAGE" \
  > "$WORK/IMAGE" 2>/dev/null || echo "unknown $APP_IMAGE" > "$WORK/IMAGE"

# ---- 4. checksums, then publish the snapshot atomically -----------------
(cd "$WORK" && find . -type f ! -name MANIFEST.sha256 -print0 | sort -z |
   xargs -0 -r sha256sum > MANIFEST.sha256)
mv "$WORK" "$FINAL"
ln -sfn "snapshots/$STAMP" "$LATEST.tmp"
mv -T "$LATEST.tmp" "$LATEST"
log "snapshot $STAMP written"

# ---- 5. keep KEEP_DAILY days and one per week for KEEP_WEEKLY weeks ----
mapfile -t snaps < <(find "$SNAP_DIR" -mindepth 1 -maxdepth 1 -type d \
                       ! -name '*.partial' -printf '%f\n' | sort -r)
declare -A week_kept=()
now_s=$(date +%s)
for s in "${snaps[@]}"; do
  day="${s%%_*}"
  age_days=$(( (now_s - $(date -d "$day" +%s)) / 86400 ))
  week="$(date -d "$day" +%G-%V)"
  if [ "$s" = "$STAMP" ] || [ "$age_days" -lt "$KEEP_DAILY" ]; then
    week_kept[$week]=1
    continue
  fi
  if [ -z "${week_kept[$week]:-}" ] && [ "$age_days" -lt $(( KEEP_WEEKLY * 7 )) ]; then
    week_kept[$week]=1
    continue
  fi
  rm -rf "${SNAP_DIR:?}/$s"
  log "pruned $s"
done
# Leftovers of runs that died before publishing.
find "$SNAP_DIR" -mindepth 1 -maxdepth 1 -type d -name '*.partial' ! -name "$STAMP.partial" \
  -mmin +60 -exec rm -rf {} + 2>/dev/null || true

# ---- 6. off the host ------------------------------------------------------
if [ -n "$BACKUP_REMOTE" ]; then
  # -H keeps the hard links, so the remote holds the same history at the
  # same cost. --delete-after mirrors the pruning above.
  rsync -aH --delete-after "$SNAP_DIR/" "$BACKUP_REMOTE/snapshots/"
  log "copied to $BACKUP_REMOTE"
else
  log "BACKUP_REMOTE is not set: this backup is on the same host as the data"
fi

date -u +%Y-%m-%dT%H:%M:%SZ > "$BACKUP_ROOT/last_success"
log "done"
