#!/usr/bin/env bash
#
# Show who owns which storage directory.
#
#   sudo deploy/scripts/list_users.sh
#   sudo deploy/scripts/list_users.sh --link      # also refresh by-name/
#
# Directories are named after the Keycloak `sub` so that a user can
# rename themselves without losing their projects. The price is that
# `ls /srv/omicsapp/users` tells you nothing, which matters the moment
# you are chasing a quota or restoring one person's backup.
#
# This asks Keycloak rather than keeping a mapping file, because a
# mapping file would be wrong the first time somebody renames themselves
# and nobody would notice until it mattered.
#
# --link maintains /srv/omicsapp/by-name/<username> -> ../users/<sub>.
# Symlinks are rebuilt from scratch each run, so a rename is picked up
# by running it again. Nothing reads them; they exist for you.
#
#   OMICSAPP_USERS_ROOT   default /srv/omicsapp/users
#   KC_URL                default http://127.0.0.1:8180/auth
#   KC_REALM              default omicsapp
#   KC_ADMIN_USER         default admin

set -euo pipefail

USERS_ROOT="${OMICSAPP_USERS_ROOT:-/srv/omicsapp/users}"
BY_NAME="${OMICSAPP_BY_NAME:-$(dirname "$USERS_ROOT")/by-name}"
KC_URL="${KC_URL:-http://127.0.0.1:8180/auth}"
KC_REALM="${KC_REALM:-omicsapp}"
KC_ADMIN_USER="${KC_ADMIN_USER:-admin}"

LINK=0
case "${1:-}" in
    --link) LINK=1 ;;
    "")     ;;
    *)      echo "usage: $0 [--link]" >&2; exit 2 ;;
esac

read -rsp "Keycloak admin password for '${KC_ADMIN_USER}': " KC_ADMIN_PASSWORD; echo

KC_ADMIN_PASSWORD="$KC_ADMIN_PASSWORD" \
KC_URL="$KC_URL" KC_REALM="$KC_REALM" KC_ADMIN_USER="$KC_ADMIN_USER" \
USERS_ROOT="$USERS_ROOT" BY_NAME="$BY_NAME" LINK="$LINK" \
python3 <<'PY'
import json, os, subprocess, sys, urllib.error, urllib.parse, urllib.request

KC_URL = os.environ["KC_URL"].rstrip("/")
REALM  = os.environ["KC_REALM"]
ROOT   = os.environ["USERS_ROOT"]
BY_NAME = os.environ["BY_NAME"]
LINK   = os.environ["LINK"] == "1"


def call(path, token=None, form=None):
    data = urllib.parse.urlencode(form).encode() if form else None
    headers = {"Content-Type": "application/x-www-form-urlencoded"} if form else {}
    if token:
        headers["Authorization"] = "Bearer " + token
    req = urllib.request.Request(KC_URL + path, data=data, headers=headers)
    with urllib.request.urlopen(req) as resp:
        raw = resp.read()
        return json.loads(raw) if raw else None


try:
    tok = call("/realms/master/protocol/openid-connect/token", form={
        "client_id": "admin-cli", "grant_type": "password",
        "username": os.environ["KC_ADMIN_USER"],
        "password": os.environ["KC_ADMIN_PASSWORD"]})
except urllib.error.HTTPError as e:
    sys.exit("error: could not authenticate to Keycloak (HTTP %d)" % e.code)
except urllib.error.URLError as e:
    sys.exit("error: could not reach Keycloak at %s (%s)" % (KC_URL, e))
token = tok["access_token"]

users = call("/admin/realms/%s/users?max=1000" % REALM, token=token) or []


def size_of(path):
    if not os.path.isdir(path):
        return "MISSING"
    try:
        out = subprocess.run(["du", "-sh", path], capture_output=True, text=True,
                             timeout=120)
        return out.stdout.split("\t")[0].strip() or "?"
    except Exception:
        return "?"


rows = []
for u in sorted(users, key=lambda x: x.get("username", "")):
    sub = u["id"]
    groups = call("/admin/realms/%s/users/%s/groups" % (REALM, sub), token=token) or []
    rows.append((u.get("username", ""), u.get("email", "") or "",
                 ",".join(g["name"] for g in groups) or "-",
                 sub, size_of(os.path.join(ROOT, sub)),
                 "" if u.get("enabled", True) else " (disabled)"))

w = max([len(r[0]) for r in rows] + [8])
print("%-*s  %-26s  %-12s  %-36s  %s" % (w, "USERNAME", "EMAIL", "GROUPS", "STORAGE ID", "SIZE"))
for username, email, groups, sub, size, flag in rows:
    print("%-*s  %-26s  %-12s  %-36s  %s%s" % (w, username, email, groups, sub, size, flag))

# A directory with no account is what a deleted user leaves behind. It
# is never removed automatically: it holds the only copy of that
# person's work, and deciding it is worthless is not this script's call.
known = {u["id"] for u in users}
if os.path.isdir(ROOT):
    orphans = [d for d in sorted(os.listdir(ROOT))
               if os.path.isdir(os.path.join(ROOT, d)) and d not in known]
    if orphans:
        print("\n%d directory(ies) with no account in Keycloak:" % len(orphans))
        for d in orphans:
            print("  %s  %s" % (d, size_of(os.path.join(ROOT, d))))

if LINK:
    os.makedirs(BY_NAME, exist_ok=True)
    for name in os.listdir(BY_NAME):
        p = os.path.join(BY_NAME, name)
        if os.path.islink(p):
            os.unlink(p)
    n = 0
    for username, _, _, sub, _, _ in rows:
        if not username or "/" in username:
            continue
        target = os.path.join(os.path.relpath(ROOT, BY_NAME), sub)
        os.symlink(target, os.path.join(BY_NAME, username))
        n += 1
    print("\n%d symlink(s) refreshed in %s" % (n, BY_NAME))
PY
