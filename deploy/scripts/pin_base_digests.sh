#!/usr/bin/env bash
#
# Pin the third-party images this deployment builds on to their digests.
#
#   deploy/scripts/pin_base_digests.sh            # resolve, record, rewrite
#   deploy/scripts/pin_base_digests.sh --check    # report drift, change nothing
#
# A tag is a name someone else can move. `bioconductor_docker:RELEASE_3_20`
# or `postgres:16` today and the same tag next month can be different
# images, so a rebuild -- or a restore onto a new host -- quietly gets
# different software from the one that was tested. A digest names the
# bytes. This script asks each registry what the tag points at now and
# writes the answer in two places:
#
#   * deploy/docker/base-digests.lock -- "<image:tag> <sha256:...>" per
#     line, the record test-deploy-contract.R checks the files against;
#   * the references themselves: `ARG BASE_IMAGE=` in the Dockerfile and
#     the `image:` lines of the Keycloak compose template, as
#     `image:tag@sha256:...` (Docker uses the digest; the tag stays for
#     the reader).
#
# Run it with network access to the registries, when you choose to take
# an update -- a security release of Postgres, a rebuilt base image --
# then rebuild and test as for any change. An image whose registry
# cannot be reached keeps its current reference and is reported; the
# contract test warns about every reference without a recorded digest.
#
# Needs only curl: it reads the registries' HTTP API directly, so it
# runs where Docker does not.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPLOY_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOCKERFILE="$DEPLOY_DIR/docker/Dockerfile"
COMPOSE="$DEPLOY_DIR/keycloak/docker-compose.yml.template"
LOCK="$DEPLOY_DIR/docker/base-digests.lock"
CURL="${CURL:-curl}"

CHECK=0
case "${1:-}" in
    --check) CHECK=1 ;;
    "") ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "pin_base_digests.sh: unknown argument: $1" >&2; exit 2 ;;
esac

ACCEPT="application/vnd.oci.image.index.v1+json, application/vnd.docker.distribution.manifest.list.v2+json, application/vnd.oci.image.manifest.v1+json, application/vnd.docker.distribution.manifest.v2+json"

# The digest a registry reports for image:tag, or nothing.
resolve() {
    local ref="$1" name tag registry repo api hdrs auth realm service scope token
    name="${ref%:*}"; tag="${ref##*:}"
    case "${name%%/*}" in
        *.*|*:*|localhost) registry="${name%%/*}"; repo="${name#*/}" ;;
        *) registry="docker.io"; repo="$name" ;;
    esac
    if [ "$registry" = "docker.io" ]; then
        api="registry-1.docker.io"
        case "$repo" in */*) ;; *) repo="library/$repo" ;; esac
    else
        api="$registry"
    fi
    hdrs="$("$CURL" -sSI -m 30 -H "Accept: $ACCEPT" "https://$api/v2/$repo/manifests/$tag" 2>/dev/null || true)"
    if printf '%s' "$hdrs" | grep -qiE '^HTTP/[0-9.]+ 401'; then
        # Anonymous pull token, from wherever the registry says to ask.
        auth="$(printf '%s' "$hdrs" | tr -d '\r' | grep -i '^www-authenticate:' | head -n 1)"
        realm="$(printf '%s' "$auth" | sed -nE 's/.*realm="([^"]+)".*/\1/p')"
        service="$(printf '%s' "$auth" | sed -nE 's/.*service="([^"]+)".*/\1/p')"
        scope="$(printf '%s' "$auth" | sed -nE 's/.*scope="([^"]+)".*/\1/p')"
        [ -n "$scope" ] || scope="repository:$repo:pull"
        [ -n "$realm" ] || return 0
        token="$("$CURL" -sS -m 30 -G "$realm" --data-urlencode "service=$service" \
                   --data-urlencode "scope=$scope" 2>/dev/null |
                 sed -nE 's/.*"(access_)?token"[[:space:]]*:[[:space:]]*"([^"]+)".*/\2/p' | head -n 1)"
        [ -n "$token" ] || return 0
        hdrs="$("$CURL" -sSI -m 30 -H "Accept: $ACCEPT" -H "Authorization: Bearer $token" \
                  "https://$api/v2/$repo/manifests/$tag" 2>/dev/null || true)"
    fi
    printf '%s' "$hdrs" | tr -d '\r' | sed -nE 's/^[Dd]ocker-[Cc]ontent-[Dd]igest:[[:space:]]*(sha256:[0-9a-f]{64}).*$/\1/p' | head -n 1
}

# Every reference: "file|image:tag" (any digest already there stripped).
references() {
    sed -nE 's/^ARG BASE_IMAGE=([^@[:space:]]+)(@sha256:[0-9a-f]+)?[[:space:]]*$/\1/p' "$DOCKERFILE" |
        sed "s|^|$DOCKERFILE\||"
    sed -nE 's/^[[:space:]]+image:[[:space:]]*([^@[:space:]]+)(@sha256:[0-9a-f]+)?[[:space:]]*$/\1/p' "$COMPOSE" |
        sed "s|^|$COMPOSE\||"
}

new_lock=""
drift=0
unresolved=0
while IFS='|' read -r file ref; do
    [ -n "$ref" ] || continue
    recorded="$( [ -f "$LOCK" ] && awk -v r="$ref" '$1 == r { print $2 }' "$LOCK" || true)"
    digest="$(resolve "$ref")"
    if [ -z "$digest" ]; then
        echo "  $ref: registry not reachable or tag unknown; kept as it is${recorded:+ (recorded $recorded)}" >&2
        unresolved=$((unresolved + 1))
        [ -n "$recorded" ] && new_lock+="$ref $recorded"$'\n'
        continue
    fi
    if [ "$digest" != "$recorded" ]; then
        drift=$((drift + 1))
        echo "  $ref: ${recorded:-unpinned} -> $digest"
    else
        echo "  $ref: $digest (unchanged)"
    fi
    new_lock+="$ref $digest"$'\n'
    if [ "$CHECK" = 0 ]; then
        esc_ref="$(printf '%s' "$ref" | sed 's/[.[\*^$/]/\\&/g')"
        sed -i -E \
            -e "s/^(ARG BASE_IMAGE=)$esc_ref(@sha256:[0-9a-f]+)?[[:space:]]*$/\1$esc_ref@$digest/" \
            -e "s/^([[:space:]]+image:[[:space:]]*)$esc_ref(@sha256:[0-9a-f]+)?[[:space:]]*$/\1$esc_ref@$digest/" \
            "$file"
    fi
done < <(references)

if [ "$CHECK" = 1 ]; then
    echo "$drift reference(s) differ from what their registry now serves; $unresolved could not be resolved."
    exit 0
fi

{
    echo "# Written by deploy/scripts/pin_base_digests.sh -- do not edit by hand."
    echo "# <image:tag> <digest>, as the registry served it on $(date -u +%Y-%m-%d)."
    printf '%s' "$new_lock"
} > "$LOCK"
echo "Wrote $LOCK. Rebuild, test, and commit the Dockerfile, compose template and lock together."
