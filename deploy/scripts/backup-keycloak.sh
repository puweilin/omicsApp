#!/bin/bash
#
# Nightly dump of the Keycloak account database.
#
# This was one crontab line:
#
#   docker exec keycloak-db pg_dump ... | gzip > /backup/.../keycloak-$(date +%u).sql.gz
#
# which fails silently in two ways at once. Cron runs it under /bin/sh, so
# there is no `pipefail` and the pipeline's exit status is gzip's -- gzip
# succeeds at compressing nothing. And the redirection truncates the target
# *before* the pipeline runs, so a night when pg_dump cannot start does not
# merely fail to produce a backup, it destroys the one from a week ago.
# What remains is a 20-byte file that `gzip -t` calls valid.
#
# So: fail on any error, write somewhere else first, prove the dump is real,
# and only then take the place of the previous one.
#
# Install:
#   sudo install -m 0755 backup-keycloak.sh /usr/local/sbin/
#   sudo mkdir -p /backup/omicsapp && sudo chmod 0700 /backup/omicsapp

set -Eeuo pipefail

CONTAINER="${KEYCLOAK_DB_CONTAINER:-keycloak-db}"
DB_USER="${KEYCLOAK_DB_USER:-keycloak}"
DB_NAME="${KEYCLOAK_DB_NAME:-keycloak}"
DEST_DIR="${BACKUP_DIR:-/backup/omicsapp}"
# Named by day of week, so seven days rotate with no cleanup job to forget.
DEST="${DEST_DIR}/keycloak-$(date +%u).sql.gz"

# A dump smaller than this is not a database. The real one is megabytes;
# the empty gzip that used to overwrite good backups is 20 bytes. This is a
# coarse screen, not proof of restorability -- only a restore is that, and
# the drill is a separate, mandatory step.
MIN_BYTES="${BACKUP_MIN_BYTES:-4096}"

# %F/%T rather than -Is: the latter is a GNU extension, and this script is
# exercised from the test suite on machines whose date(1) is BSD.
log() { printf '%s backup-keycloak: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*" >&2; }

fail() {
  log "FAILED: $*"
  # Left for the caller to route: cron mail, a systemd OnFailure unit, or
  # whatever the host actually watches. Failing quietly is how a week of
  # missing backups goes unnoticed.
  if [ -n "${BACKUP_ALERT_CMD:-}" ]; then
    "${BACKUP_ALERT_CMD}" "omicsapp keycloak backup failed: $*" || true
  fi
  exit 1
}

trap 'fail "unexpected error on line ${LINENO}"' ERR

mkdir -p "${DEST_DIR}"
chmod 0700 "${DEST_DIR}"

tmp="$(mktemp "${DEST}.XXXXXX")"
chmod 0600 "${tmp}"
trap 'rm -f "${tmp}"' EXIT

# `pipefail` is why this is a script and not a crontab line: without it the
# status below is gzip's, and gzip is perfectly happy compressing nothing.
if ! docker exec "${CONTAINER}" pg_dump -U "${DB_USER}" "${DB_NAME}" | gzip > "${tmp}"; then
  fail "pg_dump or gzip returned non-zero; the previous ${DEST} is untouched"
fi

gzip -t "${tmp}" || fail "the dump did not survive compression"

size="$(gunzip -c "${tmp}" | wc -c | tr -d ' ')"
if [ "${size}" -lt "${MIN_BYTES}" ]; then
  fail "dump is ${size} bytes, under the ${MIN_BYTES}-byte floor; the previous ${DEST} is untouched"
fi

# The dump has to contain the schema, not just be large. A Postgres error
# page compresses to something plausible too.
gunzip -c "${tmp}" | grep -q "PostgreSQL database dump" \
  || fail "the dump does not look like pg_dump output; the previous ${DEST} is untouched"

# Only now. Rename is atomic within a filesystem, so a reader either sees
# last week's backup or this one, never a partial file.
mv "${tmp}" "${DEST}"
trap - EXIT
log "wrote ${DEST} (${size} bytes uncompressed)"
