#!/usr/bin/env bash
#
# Provision a user's private project directory and print the entry to
# paste into application.yml.
#
#   sudo deploy/scripts/add_user.sh alice
#
# The directory has to exist, with the right owner, *before* the user
# first logs in. ShinyProxy bind-mounts it into the container, and
# Docker creates a missing bind-mount source owned by root -- the
# container runs as uid 1001 and would silently fail to write, leaving
# the user unable to save anything.

set -euo pipefail

USERS_ROOT="${OMICSAPP_USERS_ROOT:-/srv/omicsapp/users}"
CONTAINER_UID=1001
CONTAINER_GID=1001

usage() {
    echo "usage: $0 <username>" >&2
    echo "  creates ${USERS_ROOT}/<username> owned by ${CONTAINER_UID}:${CONTAINER_GID}" >&2
    exit 2
}

[ $# -eq 1 ] || usage
USERNAME="$1"

# The name becomes a directory under a shared root, so refuse anything
# that could escape it.
case "$USERNAME" in
    ""|*/*|*..*|.*) echo "error: invalid username '$USERNAME'" >&2; exit 1 ;;
esac

TARGET="${USERS_ROOT}/${USERNAME}"

if [ -e "$TARGET" ]; then
    echo "note: ${TARGET} already exists; leaving it untouched."
else
    install -d -o "$CONTAINER_UID" -g "$CONTAINER_GID" -m 700 "$TARGET"
    echo "created ${TARGET} (owner ${CONTAINER_UID}:${CONTAINER_GID}, mode 700)"
fi

echo
echo "Add this to the proxy.users list in application.yml:"
echo
if command -v htpasswd >/dev/null 2>&1; then
    echo "  Generate the hash with:"
    echo "    htpasswd -bnBC 10 \"\" 'THE-PASSWORD' | tr -d ':\\n'"
else
    echo "  htpasswd not found (apt-get install apache2-utils) --"
    echo "  generate the bcrypt hash on another machine."
fi
cat <<YAML

    - name: ${USERNAME}
      password: "{bcrypt}REPLACE_WITH_HASH"
      groups: [lab]

YAML
echo "Then restart ShinyProxy:  sudo systemctl restart shinyproxy"
echo "Note: a restart drops sessions that are currently running, so do"
echo "it when nobody is mid-analysis."
