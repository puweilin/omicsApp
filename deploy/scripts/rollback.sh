#!/usr/bin/env bash
#
# Switch the image ShinyProxy starts for users, and restart it.
#
#   sudo deploy/scripts/rollback.sh --list               # what is available
#   sudo deploy/scripts/rollback.sh omicsapp:0.0.0.9000-1a2b3c4d5e6f
#   sudo deploy/scripts/rollback.sh --yes <tag>          # no prompt (scripts)
#
# Named for the day it matters, and used for every image change: putting
# a new build live is the same operation as going back to an old one.
# It
#
#   1. checks the image exists on this host (an image that is not there
#      would turn the restart into an outage),
#   2. copies application.yml aside, then rewrites its one
#      `container-image:` line,
#   3. appends "<time> <old> -> <new>" to image-history beside it, so the
#      previous image is always one lookup away,
#   4. restarts ShinyProxy and waits for it to answer.
#
# Restarting ShinyProxy ends every running session -- the app autosaves,
# so users reload and restore -- hence the prompt.
#
# Projects need no rollback of their own: files written by a newer
# release stay readable by an older one (the .omp signature trailer is
# ignored, and the format version is only raised together with
# schema_version when an older reader could not cope). Rolling *forward*
# past the release that introduced signing, files the older release
# saved meanwhile are unsigned; deploy/README.md ("Rolling back") says
# how to let the store adopt them.
#
# Environment (for testing, or a non-standard layout): SP_CONF, DOCKER,
# SYSTEMCTL, SP_URL, WAIT_SECONDS, IMAGE_REPO.

set -euo pipefail

SP_CONF="${SP_CONF:-/etc/shinyproxy/application.yml}"
DOCKER="${DOCKER:-docker}"
SYSTEMCTL="${SYSTEMCTL:-systemctl}"
SP_URL="${SP_URL:-http://127.0.0.1:8081/}"
WAIT_SECONDS="${WAIT_SECONDS:-90}"
IMAGE_REPO="${IMAGE_REPO:-omicsapp}"
HISTORY="$(dirname "$SP_CONF")/image-history"

die() { echo "rollback.sh: $*" >&2; exit 1; }

current_image() {
    sed -nE 's/^[[:space:]]*container-image:[[:space:]]*"?([^"#[:space:]]+)"?.*$/\1/p' "$SP_CONF" | head -n 1
}

YES=0
TARGET=""
while [ $# -gt 0 ]; do
    case "$1" in
        -y|--yes) YES=1; shift ;;
        --list)
            [ -r "$SP_CONF" ] && echo "Running now: $(current_image)"
            echo "Images on this host (newest first):"
            "$DOCKER" image ls "$IMAGE_REPO" --format '  {{.Repository}}:{{.Tag}}\t{{.CreatedSince}}\t{{.ID}}'
            if [ -r "$HISTORY" ]; then
                echo "Recent switches ($HISTORY):"
                tail -n 5 "$HISTORY" | sed 's/^/  /'
            fi
            exit 0 ;;
        -h|--help) sed -n '2,36p' "$0"; exit 0 ;;
        -*) die "unknown option: $1" ;;
        *) [ -z "$TARGET" ] || die "one image at a time"; TARGET="$1"; shift ;;
    esac
done

[ -n "$TARGET" ] || die "which image? See: $0 --list"
case "$TARGET" in
    *:latest|latest) die "'latest' is not something you can roll back to; give an immutable tag ($0 --list)" ;;
    *[!A-Za-z0-9._:/@-]*) die "not an image reference: '$TARGET'" ;;
    # A tag (name:tag) or a digest (name@sha256:...) -- both contain a colon.
    *:*) ;;
    *) die "give the full reference with its tag, e.g. ${IMAGE_REPO}:<version>-<commit>" ;;
esac
[ -f "$SP_CONF" ] || die "no ShinyProxy configuration at $SP_CONF"
[ "$(grep -cE '^[[:space:]]*container-image:' "$SP_CONF")" = 1 ] \
    || die "expected exactly one container-image line in $SP_CONF"

if ! "$DOCKER" image inspect "$TARGET" >/dev/null 2>&1; then
    die "the image $TARGET is not on this host. Load it first, e.g.:
  gunzip -c /srv/backup/<archive>.tar.gz | sudo docker load
and see what is here with: $0 --list"
fi

CURRENT="$(current_image)"
if [ "$CURRENT" = "$TARGET" ]; then
    echo "ShinyProxy already runs $TARGET; nothing to do."
    exit 0
fi

echo "ShinyProxy runs:  ${CURRENT:-<none>}"
echo "Switching to:     $TARGET"
echo "Restarting ShinyProxy ends every running session (the app autosaves; users reload)."
if [ "$YES" != 1 ]; then
    [ -t 0 ] || die "not a terminal; pass --yes to confirm"
    read -r -p "Continue? [y/N] " answer
    case "$answer" in y|Y|yes|YES) ;; *) echo "Nothing changed."; exit 1 ;; esac
fi

STAMP="$(date +%Y%m%d-%H%M%S)"
cp -p "$SP_CONF" "$SP_CONF.bak-$STAMP"
tmp="$(mktemp "$SP_CONF.XXXXXX")"
trap 'rm -f "$tmp"' EXIT
# Only the value changes; indentation, comments and the rest stay.
sed -E "s|^([[:space:]]*container-image:[[:space:]]*)\"?[^\"#[:space:]]+\"?|\1$TARGET|" "$SP_CONF" > "$tmp"
[ "$(sed -nE 's/^[[:space:]]*container-image:[[:space:]]*([^#[:space:]]+).*$/\1/p' "$tmp")" = "$TARGET" ] \
    || die "could not rewrite container-image; $SP_CONF is unchanged"
# Same owner and mode as before: ShinyProxy runs as `shinyproxy` and must
# still be able to read a file that holds a secret (README, step 6).
chmod --reference="$SP_CONF" "$tmp" 2>/dev/null || chmod 600 "$tmp"
chown --reference="$SP_CONF" "$tmp" 2>/dev/null || true
mv "$tmp" "$SP_CONF"
trap - EXIT
printf '%s %s -> %s\n' "$(date '+%F %T')" "${CURRENT:-<none>}" "$TARGET" >> "$HISTORY"

"$SYSTEMCTL" restart shinyproxy
for _ in $(seq 1 "$WAIT_SECONDS"); do
    if curl -fsS -o /dev/null -m 5 "$SP_URL" 2>/dev/null || \
       curl -sS -o /dev/null -m 5 -w '%{http_code}' "$SP_URL" 2>/dev/null | grep -qE '^(2|3|401)'; then
        echo "ShinyProxy is up on $TARGET."
        echo "The next login starts the new image; sessions open now were ended."
        echo "To undo: sudo $0 $CURRENT"
        exit 0
    fi
    sleep 1
done
die "ShinyProxy did not answer on $SP_URL within ${WAIT_SECONDS}s. Look at:
  sudo journalctl -u shinyproxy -n 100
and to undo: sudo $0 --yes $CURRENT   (the old file is $SP_CONF.bak-$STAMP)"
