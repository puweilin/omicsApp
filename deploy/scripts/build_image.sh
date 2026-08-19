#!/usr/bin/env bash
#
# Build the omicsApp image from the repository root.
#
#   deploy/scripts/build_image.sh [tag]
#
# The build context must be the repository root: the Dockerfile COPYs
# both packages/omicsCore and packages/omicsApp. Running `docker build`
# from deploy/docker/ instead will fail on those COPY lines.

set -euo pipefail

TAG="${1:-omicsapp:1.0}"

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

if [ -n "$(git status --porcelain packages/ 2>/dev/null)" ]; then
    echo "warning: packages/ has uncommitted changes; the image will not"
    echo "         correspond to any commit." >&2
fi

echo "Building ${TAG} from ${REPO_ROOT}"
echo "Container uid: ${APP_UID}$([ "$APP_UID" = "$(id -u)" ] && echo ' (matches yours)')"
echo "First build takes 30-60 minutes (the Bioconductor stack compiles);"
echo "later builds reuse everything up to the package COPY."
echo

docker build \
  --build-arg "APP_UID=${APP_UID}" \
  -t "$TAG" -f deploy/docker/Dockerfile .

echo
echo "Built ${TAG}"
docker image inspect "$TAG" --format '  size: {{.Size}} bytes' 2>/dev/null || true
echo
echo "Smoke-test it without ShinyProxy:"
echo "  docker run --rm -p 3838:3838 -e OMICSAPP_DATA_DIR=/tmp/data ${TAG}"
echo "  then open http://localhost:3838"
