#!/usr/bin/env bash
#
# Download shellcheck, hadolint and actionlint at fixed versions, check
# each against the SHA-256 recorded here, and put them in $1 (default
# ./.lint-tools/bin), which is printed so a workflow can add it to PATH.
#
#   .github/scripts/install-lint-tools.sh "$RUNNER_TEMP/lint-tools"
#
# Fixed binaries rather than whatever the runner image or a marketplace
# action carries: a new linter release then cannot fail a build that
# changed nothing, and a tampered download fails the checksum instead of
# running in CI. To move a version, change it here together with its
# checksum (the release's own checksum file where it publishes one).

set -euo pipefail

DEST="${1:-.lint-tools/bin}"
mkdir -p "$DEST"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

SHELLCHECK_VERSION=v0.11.0
SHELLCHECK_SHA256=8c3be12b05d5c177a04c29e3c78ce89ac86f1595681cab149b65b97c4e227198
HADOLINT_VERSION=v2.15.1
HADOLINT_SHA256=c7187db94eeeeca956519a6af171adc31453941a1e777961f6e680f697c8c507
ACTIONLINT_VERSION=1.7.12
ACTIONLINT_SHA256=8aca8db96f1b94770f1b0d72b6dddcb1ebb8123cb3712530b08cc387b349a3d8

# fetch <url> <sha256> <file>
fetch() {
    curl -fsSL --retry 3 -o "$WORK/$3" "$1"
    echo "$2  $WORK/$3" | sha256sum -c --quiet - \
        || { echo "install-lint-tools.sh: checksum mismatch for $1" >&2; exit 1; }
}

fetch "https://github.com/koalaman/shellcheck/releases/download/${SHELLCHECK_VERSION}/shellcheck-${SHELLCHECK_VERSION}.linux.x86_64.tar.xz" \
    "$SHELLCHECK_SHA256" shellcheck.tar.xz
tar -xJf "$WORK/shellcheck.tar.xz" -C "$WORK"
install -m 0755 "$WORK/shellcheck-${SHELLCHECK_VERSION}/shellcheck" "$DEST/shellcheck"

fetch "https://github.com/hadolint/hadolint/releases/download/${HADOLINT_VERSION}/hadolint-linux-x86_64" \
    "$HADOLINT_SHA256" hadolint
install -m 0755 "$WORK/hadolint" "$DEST/hadolint"

fetch "https://github.com/rhysd/actionlint/releases/download/v${ACTIONLINT_VERSION}/actionlint_${ACTIONLINT_VERSION}_linux_amd64.tar.gz" \
    "$ACTIONLINT_SHA256" actionlint.tar.gz
tar -xzf "$WORK/actionlint.tar.gz" -C "$WORK" actionlint
install -m 0755 "$WORK/actionlint" "$DEST/actionlint"

cd "$DEST" && pwd
