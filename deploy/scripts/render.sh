#!/usr/bin/env bash
#
# Write the deployment files that carry the server's address.
#
#   deploy/scripts/render.sh                 # address from deploy/host.env
#   deploy/scripts/render.sh --host omics.example.org
#   deploy/scripts/render.sh --out /tmp/x    # for looking, not deploying
#
# Four files name the address, and they have to agree: nginx's
# server_name and the certificate it expects, Keycloak's KC_HOSTNAME,
# the realm's redirect URIs, ShinyProxy's OIDC endpoints. Each is kept
# as a .template with @OMICSAPP_HOST@ where the address goes, and this
# script writes the real file beside it. The rendered files are
# gitignored and excluded from the README's rsync, so they live on the
# server and survive a sync; the templates are what git tracks.
#
# @OMICSAPP_HOST_SAN@ is derived: `IP:<address>` for an IPv4 address,
# `DNS:<name>` otherwise, because Java validates a certificate against
# its subjectAltName list and the entry type has to match.
#
# Nothing here reads a secret. application.yml still needs its
# client-secret filled in on the server (keycloak/README.md).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPLOY_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

die() { echo "render.sh: $*" >&2; exit 1; }

HOST="${OMICSAPP_HOST:-}"
OUT="$DEPLOY_DIR"
while [ $# -gt 0 ]; do
    case "$1" in
        --host) [ $# -ge 2 ] || die "--host needs a value"; HOST="$2"; shift 2 ;;
        --out)  [ $# -ge 2 ] || die "--out needs a directory"; OUT="$2"; shift 2 ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

if [ -z "$HOST" ]; then
    ENV_FILE="$DEPLOY_DIR/host.env"
    [ -f "$ENV_FILE" ] || die "no $ENV_FILE. Run: cp deploy/host.env.template deploy/host.env, then set OMICSAPP_HOST in it."
    # Only the one assignment is read; the file is not sourced, so a
    # stray line in it cannot run anything.
    HOST="$(sed -nE 's/^[[:space:]]*OMICSAPP_HOST=[[:space:]]*"?([^"#[:space:]]*)"?.*$/\1/p' "$ENV_FILE" | tail -n 1)"
    [ -n "$HOST" ] && [ "$HOST" != "REPLACE_ME" ] || die "OMICSAPP_HOST is not set in $ENV_FILE"
fi

case "$HOST" in
    *[!A-Za-z0-9.-]*|"") die "OMICSAPP_HOST must be an IP address or a DNS name, without scheme, port or path: '$HOST'" ;;
esac
if [[ "$HOST" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    SAN="IP:$HOST"
else
    SAN="DNS:$HOST"
fi

# render <template relative to deploy/> <output relative to deploy/> <comment leader or ''>
render() {
    local template="$DEPLOY_DIR/$1" out="$OUT/$2" leader="$3"
    [ -f "$template" ] || die "missing template $template"
    mkdir -p "$(dirname "$out")"
    {
        if [ -n "$leader" ]; then
            printf '%s Written by deploy/scripts/render.sh from %s for %s. Edit the template, not this file.\n' \
                "$leader" "$(basename "$template")" "$HOST"
        fi
        sed -e "s|@OMICSAPP_HOST_SAN@|$SAN|g" -e "s|@OMICSAPP_HOST@|$HOST|g" "$template"
    } > "$out.tmp"
    if grep -q '@OMICSAPP_' "$out.tmp"; then
        rm -f "$out.tmp"
        die "a token in $1 has no value; this script knows @OMICSAPP_HOST@ and @OMICSAPP_HOST_SAN@"
    fi
    mv "$out.tmp" "$out"
    echo "  $out"
}

echo "Rendering for $HOST ($SAN):"
render nginx/omicsapp.conf.template            nginx/omicsapp.conf            '#'
render keycloak/docker-compose.yml.template    keycloak/docker-compose.yml    '#'
render keycloak/omicsapp-realm.json.template   keycloak/omicsapp-realm.json   ''
render shinyproxy/application.yml.template     shinyproxy/application.yml     '#'
chmod 600 "$OUT/shinyproxy/application.yml"
echo "Done. application.yml still needs client-secret (see deploy/keycloak/README.md)."
