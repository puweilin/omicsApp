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
echo "First build takes 30-60 minutes (the Bioconductor stack compiles);"
echo "later builds reuse everything up to the package COPY."
echo

docker build -t "$TAG" -f deploy/docker/Dockerfile .

echo
echo "Built ${TAG}"
docker image inspect "$TAG" --format '  size: {{.Size}} bytes' 2>/dev/null || true
echo
echo "Smoke-test it without ShinyProxy:"
echo "  docker run --rm -p 3838:3838 -e OMICSAPP_DATA_DIR=/tmp/data ${TAG}"
echo "  then open http://localhost:3838"
