#!/usr/bin/env bash
#
# Replace the `proxy.users:` list in application.yml with the entries in
# one or more generated blocks.
#
#   sudo deploy/scripts/install_users.sh omicsapp-users.yml
#   sudo deploy/scripts/install_users.sh admin.yml omicsapp-users.yml
#
# Replaces rather than appends, so running it twice is the same as
# running it once. Everything you want in the file has to be in the
# arguments -- including your own admin account. One source of truth
# beats a file that accumulates whatever each run happened to add.
#
# The original is copied aside first. Nothing is written unless the
# result parses and still holds every key the input had, because the
# failure being avoided is a config that loads with no users in it: the
# service starts, and nobody can log in.

set -euo pipefail

TARGET="${SHINYPROXY_YML:-/etc/shinyproxy/application.yml}"

[ $# -ge 1 ] || { echo "usage: $0 <users.yml> [more.yml ...]" >&2; exit 2; }
[ -f "$TARGET" ] || { echo "error: $TARGET not found" >&2; exit 1; }
for f in "$@"; do
    [ -f "$f" ] || { echo "error: $f not found" >&2; exit 1; }
done

BACKUP="${TARGET}.bak-$(date +%Y%m%d-%H%M%S)"
cp -p "$TARGET" "$BACKUP"
# A refused run must not leave a copy behind. Four rejected attempts
# should not mean four files in /etc/shinyproxy that all say the same
# thing as the one already there.
trap '[ -n "${WROTE:-}" ] || rm -f "$BACKUP"' EXIT

python3 - "$TARGET" "$BACKUP" "$@" <<'PY'
import re, sys, os, tempfile

target, backup, *sources = sys.argv[1:]
lines = open(target).read().splitlines()

# The users list lives at `  users:` under `proxy:`. Its body is
# everything indented deeper than that; the next line at the same level
# -- another proxy key, or the comment introducing one -- ends it.
start = next((i for i, l in enumerate(lines) if re.match(r"^  users:\s*$", l)), None)
if start is None:
    sys.exit("error: no '  users:' line in %s" % target)

end = len(lines)
for i in range(start + 1, len(lines)):
    if re.match(r"^  \S", lines[i]):        # 2-space key or comment
        end = i
        break
# Blank lines before that separator belong to it, not to the list.
while end > start + 1 and not lines[end - 1].strip():
    end -= 1

entries = []
for src in sources:
    for l in open(src).read().splitlines():
        if not l.strip() or l.lstrip().startswith("#"):
            continue
        entries.append(l.rstrip())

n_users = sum(1 for l in entries if re.match(r"^\s*- name:", l))
if n_users == 0:
    sys.exit("error: the input blocks contain no '- name:' entries")

# Every entry needs a hash. A REPLACE_ME reaching the live file is an
# account that exists and cannot be used, which surfaces as a login
# failure rather than as a config error.
bad = [l for l in entries if "REPLACE_ME" in l]
if bad:
    sys.exit("error: %d entr(y|ies) still say REPLACE_ME" % len(bad))

# ShinyProxy's simple backend does .password("{noop}" + user.password):
# the field is compared literally. A {bcrypt} prefix is not a hash
# there, it is a password nobody can type, and the only symptom is that
# the login fails.
enc = [l for l in entries if re.search(r'password:\s*"?\{(bcrypt|noop|pbkdf2|argon2)\}', l)]
if enc:
    sys.exit("error: %d password(s) carry an encoder prefix. Simple "
             "authentication compares the field literally, so these "
             "cannot be typed. Use the plain password." % len(enc))

n_pw = sum(1 for l in entries if re.match(r'^\s*password:', l))
if n_pw != n_users:
    sys.exit("error: %d users but %d passwords" % (n_users, n_pw))

out = lines[:start + 1] + entries + lines[end:]

# Keys outside the users list must be untouched. Counting them is a
# cheap check that the splice found the right boundaries: an `end` that
# ran too far would swallow heartbeat-rate, container-image, and the
# rest, and the file would still look plausible.
def keys(ls):
    return sorted(m.group(1) for m in
                  (re.match(r"^(\s*[a-z][a-z0-9-]*):", l) for l in ls) if m)
before = [k for k in keys(lines[:start] + lines[end:])]
after = [k for k in keys(out[:start] + out[start + 1 + len(entries):])]
if before != after:
    sys.exit("error: splice would change keys outside the users list")

# Owner as well as mode. The service runs as `shinyproxy` and the file
# is 600, so a rewrite that lands as root:root is one the service
# cannot read -- and it fails at the next restart, not at this one,
# which puts the cause and the symptom on different days.
st = os.stat(backup)
fd, tmp = tempfile.mkstemp(dir=os.path.dirname(target) or ".")
with os.fdopen(fd, "w") as fh:
    fh.write("\n".join(out) + "\n")
os.chmod(tmp, st.st_mode & 0o7777)
try:
    os.chown(tmp, st.st_uid, st.st_gid)
except PermissionError:
    pass          # not root; the mode is what mattered
os.replace(tmp, target)

print("%d user(s) written to %s" % (n_users, target))
print("previous file kept at %s" % backup)
PY

WROTE=1

echo
echo "Users now configured:"
sed -n '/^  users:/,/^  [^ ]/p' "$TARGET" | grep -E '^\s*- name:' | sed 's/^/  /'
echo
echo "Restart to apply:  sudo systemctl restart shinyproxy"
echo "A restart drops sessions that are running, so do it when nobody"
echo "is mid-analysis. If a login then fails, the previous file is at"
echo "  ${BACKUP}"
