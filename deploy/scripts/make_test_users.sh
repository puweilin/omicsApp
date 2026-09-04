#!/usr/bin/env bash
#
# Create N test accounts: private directories, random passwords, bcrypt
# hashes, and the YAML block to paste into application.yml.
#
#   sudo deploy/scripts/make_test_users.sh 10
#
# Passwords are generated here and written to a file only you can read.
# They are never printed to the terminal, so they do not end up in a
# scrollback buffer or a shell history, and nothing has to type them
# anywhere on your behalf.
#
# Two files land next to you:
#
#   omicsapp-users.yml        hashes only -- paste into application.yml
#   omicsapp-passwords.txt    mode 600 -- hand out, then delete
#
# Delete the second one once the accounts are distributed. A bcrypt
# hash is not reversible; that file is the only copy.

set -euo pipefail

COUNT="${1:-10}"
PREFIX="${USER_PREFIX:-test}"
GROUP="${USER_GROUP:-lab}"
COST="${BCRYPT_COST:-10}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ADD_USER="${SCRIPT_DIR}/add_user.sh"

case "$COUNT" in
    ''|*[!0-9]*) echo "usage: $0 <count>   (got '$COUNT')" >&2; exit 2 ;;
esac
[ "$COUNT" -ge 1 ] && [ "$COUNT" -le 200 ] || {
    echo "error: count must be 1-200" >&2; exit 2; }

command -v htpasswd >/dev/null 2>&1 || {
    echo "error: htpasswd not found. sudo apt-get install -y apache2-utils" >&2
    exit 1; }
[ -x "$ADD_USER" ] || { echo "error: $ADD_USER not found" >&2; exit 1; }

# Written where the invoking user can reach them, not in root's cwd.
OUT_DIR="${OUT_DIR:-$PWD}"
YAML="${OUT_DIR}/omicsapp-users.yml"
SECRETS="${OUT_DIR}/omicsapp-passwords.txt"

for f in "$YAML" "$SECRETS"; do
    [ -e "$f" ] && { echo "error: $f exists; move it aside first" >&2; exit 1; }
done

# 600 before anything is written to it, not after.
umask 077
: > "$SECRETS"
: > "$YAML"

printf '# Paste under `proxy.users:` in /etc/shinyproxy/application.yml\n' >> "$YAML"
printf '# Generated %s\n' "$(date -Iseconds)" >> "$YAML"
printf '%-12s  %s\n' "ACCOUNT" "PASSWORD" >> "$SECRETS"

for i in $(seq 1 "$COUNT"); do
    name=$(printf '%s%02d' "$PREFIX" "$i")

    # 12 characters from an alphabet with no quoting hazards for YAML
    # and none of the pairs that are read wrong aloud (0/O, 1/l/I).
    #
    # The bounded `head` comes first on purpose. Reading /dev/urandom
    # into `head -c 12` at the end of the pipe hands `tr` a SIGPIPE the
    # moment head has enough, which is a non-zero exit, which under
    # `set -o pipefail` ends the script -- silently, before the first
    # account is written.
    raw=$(LC_ALL=C head -c 256 /dev/urandom | base64 \
          | LC_ALL=C tr -dc 'A-HJ-NP-Za-km-z2-9')
    pw=${raw:0:12}
    [ ${#pw} -eq 12 ] || { echo "error: could not generate a password" >&2; exit 1; }

    hash=$(htpasswd -bnBC "$COST" "" "$pw" | tr -d ':\n')

    "$ADD_USER" "$name" >/dev/null

    printf '    - name: %s\n      password: "{bcrypt}%s"\n      groups: [%s]\n' \
        "$name" "$hash" "$GROUP" >> "$YAML"
    printf '%-12s  %s\n' "$name" "$pw" >> "$SECRETS"

    echo "created ${name}"
done

# Under sudo these belong to the person who ran it, not to root.
if [ -n "${SUDO_USER:-}" ]; then
    chown "$SUDO_USER" "$YAML" "$SECRETS" 2>/dev/null || true
fi

cat <<DONE

${COUNT} accounts created.

  ${YAML}
      the accounts. Install them with:
          sudo ${SCRIPT_DIR}/install_users.sh ${YAML}
      Add your own admin account to that file first -- install_users.sh
      replaces the list rather than adding to it, so whatever is not in
      the file will not be in the config.

  ${SECRETS}
      mode 600, the only copy of the passwords. Hand them out, then:
          shred -u ${SECRETS}

Restarting ShinyProxy drops sessions that are running, so do it when
nobody is mid-analysis.
DONE
