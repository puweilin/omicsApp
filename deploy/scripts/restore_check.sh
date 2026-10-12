#!/usr/bin/env bash
#
# Weekly restore drill: proves the latest backup can actually be
# restored, and complains when it cannot -- or when there is no recent
# backup at all.
#
#   sudo deploy/scripts/restore_check.sh
#
# A backup that has never been restored is a hope. This checks, without
# touching the running service:
#
#   1. a snapshot finished within MAX_AGE_HOURS (cron may have stopped,
#      the disk may be full -- backup.sh cannot report its own absence);
#   2. every file still matches the checksum written when it was taken;
#   3. the Keycloak dump restores into a throwaway Postgres container and
#      holds the accounts (user_entity) the live database holds;
#   4. a project file opens with the image's own omicsCore
#      (load_project()), so the .omp files are readable, not just present;
#   5. with BACKUP_REMOTE set, the remote copy has the same latest snapshot.
#
# Failures go through the same alert hook as backup.sh.

set -Eeuo pipefail

CONF="${OMICSAPP_BACKUP_CONF:-/etc/omicsapp/backup.env}"
# shellcheck disable=SC1090
[ -r "$CONF" ] && . "$CONF"

BACKUP_ROOT="${BACKUP_ROOT:-/backup/omicsapp}"
MAX_AGE_HOURS="${MAX_AGE_HOURS:-36}"
BACKUP_REMOTE="${BACKUP_REMOTE:-}"
KEYCLOAK_CONTAINER="${KEYCLOAK_CONTAINER:-keycloak-db}"
SP_CONF="${SP_CONF:-/etc/shinyproxy/application.yml}"
# The image ShinyProxy runs, unless told otherwise: after a rollback that
# is the one whose load_project() has to open these files.
if [ -z "${APP_IMAGE:-}" ] && [ -r "$SP_CONF" ]; then
  APP_IMAGE="$(sed -nE 's/^[[:space:]]*container-image:[[:space:]]*"?([^"#[:space:]]+)"?.*$/\1/p' "$SP_CONF" | head -n 1)"
fi
APP_IMAGE="${APP_IMAGE:-omicsapp:1.0}"
POSTGRES_IMAGE="${POSTGRES_IMAGE:-postgres:16}"
ALERT_WEBHOOK="${ALERT_WEBHOOK:-}"
ALERT_EMAIL="${ALERT_EMAIL:-}"
DOCKER="${DOCKER:-docker}"
# Set to 1 to skip the steps that need docker (used by the test suite).
SKIP_DOCKER_CHECKS="${SKIP_DOCKER_CHECKS:-0}"

failures=0
log() { printf '%s restore-check: %s\n' "$(date '+%F %T')" "$*" >&2; }
fail() {
  failures=$((failures + 1))
  local msg
  msg="omicsApp restore check on $(hostname): $*"
  log "FAIL: $*"
  command -v logger >/dev/null 2>&1 && logger -t omicsapp-backup -p user.err "$msg" || true
  if [ -n "$ALERT_WEBHOOK" ]; then
    local body
    body="$(printf '%s' "$msg" | sed 's/\\/\\\\/g; s/"/\\"/g')"
    curl -fsS -m 20 -H 'Content-Type: application/json' \
      -d "{\"text\":\"$body\"}" "$ALERT_WEBHOOK" >/dev/null 2>&1 || true
  fi
  if [ -n "$ALERT_EMAIL" ] && command -v mail >/dev/null 2>&1; then
    printf '%s\n' "$msg" | mail -s "omicsApp restore check FAILED" "$ALERT_EMAIL" || true
  fi
  if [ -n "${ALERT_LOG:-}" ]; then printf '%s\n' "$msg" >> "$ALERT_LOG"; fi
}

latest="$BACKUP_ROOT/latest"
if [ ! -d "$latest" ]; then
  fail "no snapshot at $latest"
  exit 1
fi
snap="$(readlink -f "$latest")"
log "checking $snap"

# ---- 1. recent enough -------------------------------------------------------
if [ -r "$BACKUP_ROOT/last_success" ]; then
  last="$(date -d "$(cat "$BACKUP_ROOT/last_success")" +%s)"
  age_h=$(( ($(date +%s) - last) / 3600 ))
  if [ "$age_h" -gt "$MAX_AGE_HOURS" ]; then
    fail "the last successful backup is ${age_h} h old (limit ${MAX_AGE_HOURS} h)"
  fi
else
  fail "no record of a successful backup ($BACKUP_ROOT/last_success)"
fi

# ---- 2. checksums ------------------------------------------------------------
if [ -r "$snap/MANIFEST.sha256" ]; then
  if ! (cd "$snap" && sha256sum --quiet -c MANIFEST.sha256 >/dev/null 2>&1); then
    fail "files in $snap no longer match their checksums"
  fi
else
  fail "$snap has no MANIFEST.sha256"
fi

if [ "$SKIP_DOCKER_CHECKS" != "1" ]; then
  # ---- 3. the account database restores ------------------------------------
  if [ -r "$snap/keycloak.sql.gz" ]; then
    name="omicsapp-restore-check-$$"
    cleanup() { "$DOCKER" rm -f "$name" >/dev/null 2>&1 || true; }
    trap cleanup EXIT
    "$DOCKER" run -d --rm --name "$name" -e POSTGRES_PASSWORD=check \
      -e POSTGRES_USER=keycloak -e POSTGRES_DB=keycloak "$POSTGRES_IMAGE" >/dev/null
    for _ in $(seq 1 60); do
      "$DOCKER" exec "$name" pg_isready -U keycloak >/dev/null 2>&1 && break
      sleep 1
    done
    if gzip -cd "$snap/keycloak.sql.gz" |
         "$DOCKER" exec -i "$name" psql -q -v ON_ERROR_STOP=1 -U keycloak keycloak >/dev/null; then
      restored="$("$DOCKER" exec "$name" psql -tA -U keycloak keycloak \
                    -c 'select count(*) from user_entity' 2>/dev/null || echo 0)"
      live="$("$DOCKER" exec "$KEYCLOAK_CONTAINER" psql -tA -U keycloak keycloak \
                -c 'select count(*) from user_entity' 2>/dev/null || echo "")"
      log "accounts: $restored restored, ${live:-?} live"
      if [ "${restored:-0}" -eq 0 ]; then
        fail "the restored account database has no users"
      elif [ -n "$live" ] && [ "$restored" -lt "$live" ]; then
        log "note: $((live - restored)) account(s) were added since the snapshot"
      fi
    else
      fail "the Keycloak dump did not restore"
    fi
  else
    fail "$snap has no keycloak.sql.gz"
  fi

  # ---- 4. a project opens --------------------------------------------------
  omp="$(find "$snap/users" -name '*.omp' -type f 2>/dev/null | shuf -n 1 || true)"
  if [ -n "$omp" ]; then
    # With the store's own key when it kept one (no OMICSAPP_SIGNING_KEY),
    # so the drill also proves the signature survived: a project the app
    # would refuse to open is not a restorable one.
    key_mount=()
    if [ -r "$(dirname "$omp")/.omicsapp-signing-key" ]; then
      key_mount=(-v "$(dirname "$omp")/.omicsapp-signing-key:/check.key:ro")
    fi
    # The same confinement the app gets, and no network at all: this
    # opens a file from a backup.
    if ! "$DOCKER" run --rm --network none --read-only --tmpfs /tmp \
          --cap-drop ALL --security-opt no-new-privileges --pids-limit 256 \
          -v "$omp:/check.omp:ro" "${key_mount[@]}" "$APP_IMAGE" \
          R -q -e 'k <- if (file.exists("/check.key")) readLines("/check.key", 1L) else NULL; p <- omicsCore::load_project("/check.omp", signing_key = k); stopifnot(omicsCore::is_omics_project(p))' \
          >/dev/null 2>&1; then
      fail "a backed-up project does not open: ${omp#"$snap"/}"
    fi
  else
    log "no project files to try yet"
  fi
fi

# ---- 5. the off-host copy keeps up -------------------------------------------
if [ -n "$BACKUP_REMOTE" ]; then
  want="$(basename "$snap")"
  if ! rsync --list-only "$BACKUP_REMOTE/snapshots/$want/MANIFEST.sha256" >/dev/null 2>&1; then
    fail "the remote copy at $BACKUP_REMOTE lacks the latest snapshot $want"
  fi
fi

if [ "$failures" -gt 0 ]; then
  log "$failures check(s) failed"
  exit 1
fi
log "all checks passed"
