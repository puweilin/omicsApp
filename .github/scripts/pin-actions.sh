#!/usr/bin/env bash
#
# Check, or re-pin, every third-party action in .github/workflows.
#
#   .github/scripts/pin-actions.sh            # check: each SHA is what its tag says
#   .github/scripts/pin-actions.sh --update   # rewrite each SHA from its tag
#
# Every `uses:` names a full commit SHA, with the release it came from
# in a trailing comment:
#
#   - uses: actions/checkout@11d5960a326750d5838078e36cf38b85af677262 # v4.4.0
#
# A tag can be moved to other code by whoever controls the repository; a
# commit SHA cannot. The comment is what a reader -- and Dependabot,
# which updates both together (.github/dependabot.yml) -- goes by. This
# script resolves each comment's tag with `git ls-remote` and compares,
# so a hand edit that changed one without the other is caught. --update
# writes the resolved SHA instead (to move to a new release, edit the
# comment's tag and run it). Needs only git and network access to
# github.com.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
UPDATE=0
[ "${1:-}" = "--update" ] && UPDATE=1

status=0
while IFS= read -r line; do
    file="${line%%:*}"
    rest="${line#*:}"
    # owner/repo[/path]@sha # tag
    ref="$(sed -nE 's/.*uses:[[:space:]]*([^@[:space:]]+)@([0-9a-f]{40})[[:space:]]*#[[:space:]]*(v[^[:space:]]+).*/\1 \2 \3/p' <<<"$rest")"
    if [ -z "$ref" ]; then
        echo "UNPINNED  $file: $(sed -E 's/^[[:space:]]+//' <<<"$rest")"
        status=1
        continue
    fi
    read -r action sha tag <<<"$ref"
    repo="$(cut -d/ -f1-2 <<<"$action")"
    # An annotated tag lists the tag object and, with ^{}, the commit it
    # points to; the commit is what `uses:` must name.
    resolved="$(git ls-remote "https://github.com/$repo" "refs/tags/$tag" "refs/tags/$tag^{}" \
        | awk '{print $1}' | tail -n 1)"
    if [ -z "$resolved" ]; then
        echo "NO TAG    $file: $repo $tag"
        status=1
    elif [ "$resolved" != "$sha" ]; then
        if [ "$UPDATE" = 1 ]; then
            sed -i "s|$action@$sha|$action@$resolved|g" "$file"
            echo "UPDATED   $file: $action $tag -> $resolved"
        else
            echo "MISMATCH  $file: $action@$sha but $tag is $resolved"
            status=1
        fi
    else
        echo "ok        $action $tag"
    fi
done < <(grep -HE '^[[:space:]-]*uses:[[:space:]]*[^.[:space:]]' .github/workflows/*.y*ml | grep -v 'uses: \./')

exit "$status"
