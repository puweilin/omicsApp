#!/usr/bin/env bash
#
# Provision a user's private directories and print the entry to paste
# into application.yml.
#
#   sudo deploy/scripts/add_user.sh alice
#
# The directories have to exist, with the right owner, *before* the user
# first logs in. ShinyProxy bind-mounts them into the container, and
# Docker creates a missing bind-mount source owned by root -- the
# container runs unprivileged and would silently fail to write, leaving
# the user unable to save anything and nothing on screen to say why.
#
# Two directories, because they are used differently: projects are small
# and touched on every save, archived uploads are bulk and written once.
# Set RAW_ROOT to put the second on another disk, or leave it unset to
# keep everything together.
#
#   USERS_ROOT   projects            default /srv/omicsapp/users
#   RAW_ROOT     archived uploads    unset = live under USERS_ROOT
#   CONTAINER_UID  must match the image's APP_UID build arg

set -euo pipefail

USERS_ROOT="${OMICSAPP_USERS_ROOT:-/srv/omicsapp/users}"
RAW_ROOT="${OMICSAPP_RAW_ROOT:-}"
CONTAINER_UID="${APP_UID:-1001}"
CONTAINER_GID="${APP_GID:-$CONTAINER_UID}"

usage() {
    cat >&2 <<USAGE
usage: $0 <username>

  creates ${USERS_ROOT}/<username>
$([ -n "$RAW_ROOT" ] && echo "      and ${RAW_ROOT}/<username>")
  owned by ${CONTAINER_UID}:${CONTAINER_GID}, mode 700

environment:
  OMICSAPP_USERS_ROOT   project directory root (default ${USERS_ROOT})
  OMICSAPP_RAW_ROOT     archive root on another disk (default: under USERS_ROOT)
  APP_UID               container uid (default 1001; must match the image)
USAGE
    exit 2
}

[ $# -eq 1 ] || usage
USERNAME="$1"

# The name becomes a directory under a shared root, so refuse anything
# that could escape it.
case "$USERNAME" in
    ""|*/*|*..*|.*) echo "error: invalid username '$USERNAME'" >&2; exit 1 ;;
esac

make_dir() {
    local target="$1" label="$2"
    if [ -e "$target" ]; then
        echo "note: ${target} already exists; leaving it untouched."
        return
    fi
    install -d -o "$CONTAINER_UID" -g "$CONTAINER_GID" -m 700 "$target"
    echo "created ${target}  (${label}, owner ${CONTAINER_UID}:${CONTAINER_GID}, mode 700)"
}

make_dir "${USERS_ROOT}/${USERNAME}" "projects"
if [ -n "$RAW_ROOT" ]; then
    make_dir "${RAW_ROOT}/${USERNAME}" "archived uploads"
else
    # The default, and the deployed layout: everything on one disk. The
    # app creates raw/ under the project directory on first upload, so
    # there is nothing to do here.
    echo "  archived uploads will go to ${USERS_ROOT}/${USERNAME}/raw"
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
