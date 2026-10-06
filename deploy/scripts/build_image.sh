#!/usr/bin/env bash
#
# Build the omicsApp image from the repository root, under a tag that
# names exactly what is in it.
#
#   deploy/scripts/build_image.sh                  # omicsapp:<version>-<commit>
#   deploy/scripts/build_image.sh omicsapp:1.0     # ...and an alias as well
#
# Every build is tagged <repo>:<version>-<commit>: the omicsApp package
# version and the first 12 characters of the git commit it was built
# from. That tag is never reused -- a build of a commit that already has
# one stops rather than overwrite it -- so yesterday's image is still
# there, under a name you can give deploy/scripts/rollback.sh, after
# today's goes wrong. Uncommitted changes get a `-dirty.<time>` suffix
# instead: an image that matches no commit should not look as if it did.
#
# Extra arguments are additional tags (convenient aliases such as
# omicsapp:1.0). `latest` is refused: it names whatever was built last,
# which is the one thing a rollback cannot use.
#
# Without a .git directory -- the server's copy is an rsync that leaves
# .git out -- the commit is read from REVISION at the repository root,
# which step 0 of deploy/README.md writes before copying. Neither present
# means a `-norev.<time>` suffix and a warning.
#
# The build context must be the repository root: the Dockerfile COPYs
# both packages/omicsCore and packages/omicsApp. Running `docker build`
# from deploy/docker/ instead will fail on those COPY lines.
#
# Environment: IMAGE_REPO (default omicsapp), APP_UID (default 1001),
# DOCKER (default docker).

set -euo pipefail

IMAGE_REPO="${IMAGE_REPO:-omicsapp}"
DOCKER="${DOCKER:-docker}"

# The uid the container runs as. It has to match whoever owns the
# per-user directories on the host, because the kernel compares numbers,
# not names. Leave it alone if you have root and will chown those
# directories to 1001; set APP_UID=$(id -u) if you do not, and own them
# yourself instead.
APP_UID="${APP_UID:-1001}"

# Resolve the repository root from this script's location so the build
# works regardless of where it is invoked from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

cd "$REPO_ROOT"

die() { echo "build_image.sh: $*" >&2; exit 1; }

VERSION="$(sed -nE 's/^Version:[[:space:]]*([^[:space:]]+).*$/\1/p' packages/omicsApp/DESCRIPTION | head -n 1)"
[ -n "$VERSION" ] || die "cannot read Version from packages/omicsApp/DESCRIPTION"
STAMP="$(date -u +%Y%m%d%H%M%S)"

if git -C "$REPO_ROOT" rev-parse --git-dir >/dev/null 2>&1; then
    REVISION="$(git -C "$REPO_ROOT" rev-parse HEAD)"
    if [ -n "$(git -C "$REPO_ROOT" status --porcelain packages/ deploy/docker/ 2>/dev/null)" ]; then
        echo "warning: packages/ or deploy/docker/ has uncommitted changes; the image" >&2
        echo "         will not correspond to any commit, and its tag says so." >&2
        SUFFIX="-dirty.${STAMP}"
    else
        SUFFIX=""
    fi
elif [ -f REVISION ]; then
    REVISION="$(awk 'NR == 1 { print $1 }' REVISION)"
    if grep -q dirty REVISION; then SUFFIX="-dirty.${STAMP}"; else SUFFIX=""; fi
else
    echo "warning: no .git and no REVISION file, so the commit is unknown." >&2
    echo "         Write one before copying: git rev-parse HEAD > REVISION" >&2
    REVISION="unknown"
    SUFFIX="-norev.${STAMP}"
fi
[[ "$REVISION" =~ ^([0-9a-f]{40}|unknown)$ ]] || die "REVISION does not hold a commit hash: '$REVISION'"

SHORT="${REVISION:0:12}"
[ "$REVISION" = "unknown" ] && SHORT="unknown"
IMMUTABLE="${IMAGE_REPO}:${VERSION}-${SHORT}${SUFFIX}"
[ "$REVISION" = "unknown" ] && IMMUTABLE="${IMAGE_REPO}:${VERSION}${SUFFIX}"

EXTRA=()
for tag in "$@"; do
    case "$tag" in
        *:latest|latest) die "refusing the tag '$tag': use the immutable tag ${IMMUTABLE}" ;;
        *:*) EXTRA+=("$tag") ;;
        *) die "'$tag' is not a tag (expected name:tag, e.g. omicsapp:1.0)" ;;
    esac
done

if "$DOCKER" image inspect "$IMMUTABLE" >/dev/null 2>&1; then
    echo "${IMMUTABLE} already exists; it is not rebuilt, so the tag keeps naming one image."
    echo "Nothing to do. (To rebuild anyway, remove it first: docker image rm ${IMMUTABLE})"
    for tag in "${EXTRA[@]}"; do "$DOCKER" tag "$IMMUTABLE" "$tag"; echo "  also tagged $tag"; done
    exit 0
fi

echo "Building ${IMMUTABLE} from ${REPO_ROOT}"
echo "Container uid: ${APP_UID}$([ "$APP_UID" = "$(id -u)" ] && echo ' (matches yours)')"
echo "First build takes 30-60 minutes (the Bioconductor stack compiles);"
echo "later builds reuse everything up to the package COPY."
echo

TAG_ARGS=(-t "$IMMUTABLE")
for tag in "${EXTRA[@]}"; do TAG_ARGS+=(-t "$tag"); done

"$DOCKER" build \
  --build-arg "APP_UID=${APP_UID}" \
  --label "org.opencontainers.image.version=${VERSION}" \
  --label "org.opencontainers.image.revision=${REVISION}" \
  --label "org.opencontainers.image.created=$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  "${TAG_ARGS[@]}" -f deploy/docker/Dockerfile .

echo
echo "Built ${IMMUTABLE}"
for tag in "${EXTRA[@]}"; do echo "  also tagged ${tag}"; done
"$DOCKER" image inspect "$IMMUTABLE" --format '  size: {{.Size}} bytes' 2>/dev/null || true
echo
echo "Smoke-test it without ShinyProxy:"
echo "  docker run --rm -p 3838:3838 -e OMICSAPP_DATA_DIR=/tmp/data ${IMMUTABLE}"
echo "  then open http://localhost:3838"
echo "Put it live (and the same command takes you back to any earlier tag):"
echo "  sudo deploy/scripts/rollback.sh ${IMMUTABLE}"
echo "Earlier builds stay available: docker image ls ${IMAGE_REPO}"
