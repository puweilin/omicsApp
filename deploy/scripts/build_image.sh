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

# What went into the image, recorded as a digest of the sources rather
# than as a commit.
#
# The commit is not knowable here. Deployment is an rsync, and
# rsync.exclude leaves `.git/` alone -- so on the server it holds whatever
# was cloned months ago, and `git rev-parse HEAD` answers a question about
# that clone rather than about the files just copied in. A version check
# that reads it is worse than none: it reports "clean" against history the
# copied files have nothing to do with.
#
# A digest over the sources does answer it, and it answers the same way on
# the machine the code came from.
SRC_DIGEST="$(find packages -type f \( -name '*.R' -o -name '*.Rmd' \
    -o -name 'DESCRIPTION' -o -name 'NAMESPACE' \) -print0 \
    | LC_ALL=C sort -z \
    | xargs -0 shasum -a 256 2>/dev/null \
    | shasum -a 256 | cut -c1-16)"
echo "Source digest: ${SRC_DIGEST}"
echo "  Compare with the same command where the code came from; the"
echo "  server's .git is not copied and cannot answer this."

echo "Building ${TAG} from ${REPO_ROOT}"
echo "Container uid: ${APP_UID}$([ "$APP_UID" = "$(id -u)" ] && echo ' (matches yours)')"
echo "First build takes 30-60 minutes (the Bioconductor stack compiles);"
echo "later builds reuse everything up to the package COPY."
echo

docker build \
  --build-arg "APP_UID=${APP_UID}" \
  --label "org.omicsapp.source-digest=${SRC_DIGEST}" \
  -t "$TAG" -f deploy/docker/Dockerfile .

echo
echo "Built ${TAG}"
echo "  source digest: ${SRC_DIGEST}"
echo "  ask a built image later with:"
echo "    docker image inspect ${TAG} --format '{{index .Config.Labels \"org.omicsapp.source-digest\"}}'"
docker image inspect "$TAG" --format '  size: {{.Size}} bytes' 2>/dev/null || true
echo
echo "Smoke-test it without ShinyProxy:"
echo "  docker run --rm -p 3838:3838 -e OMICSAPP_DATA_DIR=/tmp/data ${TAG}"
echo "  then open http://localhost:3838"
