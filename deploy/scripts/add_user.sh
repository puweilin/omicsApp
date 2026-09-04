#!/usr/bin/env bash
#
# Create omicsApp accounts in Keycloak and their storage on disk.
#
#   sudo deploy/scripts/add_user.sh alice@example.com bob@example.com
#
# One command per person, and application.yml is never touched. That is
# the point of moving to Keycloak: adding someone no longer means
# editing the config and restarting a service that would drop everyone
# else's running analysis.
#
# What each account gets:
#
#   username        the email address, to start with. The user may
#                   rename themselves later; nothing here depends on it.
#   group           lab -- which is what access-groups in application.yml
#                   checks. An account outside it can sign in and will
#                   see no applications.
#   password        random, temporary. Keycloak forces a change at first
#                   login, so the password written here stops working the
#                   moment the user has chosen their own. You never learn
#                   what they picked.
#   directory       /srv/omicsapp/users/<sub>, owned by the container uid
#
# The directory is named after the Keycloak `sub`, a UUID fixed when the
# account is created, because a user who renames themselves must not
# lose their projects. It is unreadable on purpose -- list_users.sh maps
# it back to a name.
#
# It has to exist, with the right owner, *before* that user first logs
# in. ShinyProxy bind-mounts it into the container, and Docker creates a
# missing bind-mount source owned by root; the container runs
# unprivileged and would silently fail to write, leaving the user unable
# to save anything and nothing on screen to say why. So the order below
# is deliberate: create the account, learn its sub, make the directory,
# and only then grant the group that lets it start an app.
#
#   OMICSAPP_USERS_ROOT   project directory root (default /srv/omicsapp/users)
#   OMICSAPP_RAW_ROOT     archive root on another disk (default: under USERS_ROOT)
#   APP_UID               container uid (default 1001; must match the image)
#   KC_URL                Keycloak base URL (default http://127.0.0.1:8180/auth)
#   KC_REALM              realm (default omicsapp)
#   KC_ADMIN_USER         Keycloak admin (default admin)

set -euo pipefail

USERS_ROOT="${OMICSAPP_USERS_ROOT:-/srv/omicsapp/users}"
RAW_ROOT="${OMICSAPP_RAW_ROOT:-}"
CONTAINER_UID="${APP_UID:-1001}"
CONTAINER_GID="${APP_GID:-$CONTAINER_UID}"

# Plain HTTP to the container's published port rather than through
# nginx: this runs on the server, and going the long way round would
# mean teaching curl to trust the self-signed certificate for no gain.
KC_URL="${KC_URL:-http://127.0.0.1:8180/auth}"
KC_REALM="${KC_REALM:-omicsapp}"
KC_ADMIN_USER="${KC_ADMIN_USER:-admin}"
KC_GROUP="${KC_GROUP:-lab}"

usage() {
    cat >&2 <<USAGE
usage: $0 <email> [email ...]

  creates each account in Keycloak realm '${KC_REALM}', adds it to the
  '${KC_GROUP}' group, and creates ${USERS_ROOT}/<sub>
  owned by ${CONTAINER_UID}:${CONTAINER_GID}, mode 700

environment:
  OMICSAPP_USERS_ROOT   project directory root (default ${USERS_ROOT})
  OMICSAPP_RAW_ROOT     archive root on another disk (default: under USERS_ROOT)
  APP_UID               container uid (default 1001; must match the image)
  KC_URL                Keycloak base URL (default ${KC_URL})
  KC_REALM              realm (default ${KC_REALM})
  KC_ADMIN_USER         Keycloak admin user (default ${KC_ADMIN_USER})
USAGE
    exit 2
}

[ $# -ge 1 ] || usage

for email in "$@"; do
    case "$email" in
        *@*.*) ;;
        *) echo "error: '$email' does not look like an email address" >&2; exit 1 ;;
    esac
    # The local part of an address may legally contain a slash. It would
    # never reach a path here -- directories are named after the sub --
    # but rejecting it costs nothing and keeps that true if this ever
    # changes.
    case "$email" in
        */*|*..*|.*) echo "error: refusing '$email'" >&2; exit 1 ;;
    esac
done

[ "$(id -u)" -eq 0 ] || {
    echo "error: run with sudo -- creating directories owned by uid ${CONTAINER_UID} needs root" >&2
    exit 1; }

command -v python3 >/dev/null || { echo "error: python3 not found" >&2; exit 1; }
command -v curl    >/dev/null || { echo "error: curl not found" >&2; exit 1; }

# Read rather than accept as an argument: an argument would sit in the
# shell history and be visible in `ps` to every other account on this
# machine -- and this one opens the whole realm.
read -rsp "Keycloak admin password for '${KC_ADMIN_USER}': " KC_ADMIN_PASSWORD; echo
[ -n "$KC_ADMIN_PASSWORD" ] || { echo "error: empty password" >&2; exit 1; }

# Written where the person who ran sudo can reach it, not in root's cwd.
OUT_DIR="${OUT_DIR:-${SUDO_USER:+/home/$SUDO_USER}}"
OUT_DIR="${OUT_DIR:-$PWD}"
SECRETS="${OUT_DIR}/omicsapp-passwords.txt"

# 600 before anything is written to it, not after.
umask 077

# The admin password goes in the environment, not in argv: /proc/*/environ
# is readable only by the process owner, while /proc/*/cmdline is world
# readable, so an argument would expose it to anyone with a shell here
# for as long as the script runs.
KC_ADMIN_PASSWORD="$KC_ADMIN_PASSWORD" \
KC_URL="$KC_URL" KC_REALM="$KC_REALM" KC_ADMIN_USER="$KC_ADMIN_USER" \
KC_GROUP="$KC_GROUP" USERS_ROOT="$USERS_ROOT" RAW_ROOT="$RAW_ROOT" \
CONTAINER_UID="$CONTAINER_UID" CONTAINER_GID="$CONTAINER_GID" \
SECRETS="$SECRETS" \
python3 - "$@" <<'PY'
import json, os, secrets, string, sys, urllib.error, urllib.parse, urllib.request

KC_URL   = os.environ["KC_URL"].rstrip("/")
REALM    = os.environ["KC_REALM"]
GROUP    = os.environ["KC_GROUP"]
UID      = int(os.environ["CONTAINER_UID"])
GID      = int(os.environ["CONTAINER_GID"])
USERS_ROOT = os.environ["USERS_ROOT"]
RAW_ROOT   = os.environ["RAW_ROOT"]
SECRETS    = os.environ["SECRETS"]

emails = sys.argv[1:]

# os.makedirs would happily create a missing root as well, owned by root
# and mode 700 -- and then every account under it would be unreachable
# for exactly the reason this script exists to prevent. A missing root
# means the storage is not mounted, which is worth stopping for.
for root in (USERS_ROOT, RAW_ROOT):
    if root and not os.path.isdir(root):
        sys.exit("error: %s does not exist. Is /srv/omicsapp mounted?" % root)


def call(method, path, token=None, body=None, form=None):
    url = KC_URL + path
    if form is not None:
        data = urllib.parse.urlencode(form).encode()
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
    elif body is not None:
        data = json.dumps(body).encode()
        headers = {"Content-Type": "application/json"}
    else:
        data, headers = None, {}
    if token:
        headers["Authorization"] = "Bearer " + token
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    with urllib.request.urlopen(req) as resp:
        raw = resp.read()
        return resp.status, resp.headers, (json.loads(raw) if raw else None)


def fail(msg, err=None):
    if isinstance(err, urllib.error.HTTPError):
        detail = err.read().decode(errors="replace")[:300]
        msg = "%s (HTTP %d) %s" % (msg, err.code, detail)
    elif err is not None:
        msg = "%s (%s)" % (msg, err)
    sys.exit("error: " + msg)


try:
    _, _, tok = call("POST", "/realms/master/protocol/openid-connect/token", form={
        "client_id": "admin-cli",
        "grant_type": "password",
        "username": os.environ["KC_ADMIN_USER"],
        "password": os.environ["KC_ADMIN_PASSWORD"],
    })
except urllib.error.HTTPError as e:
    fail("could not authenticate to Keycloak at " + KC_URL, e)
except urllib.error.URLError as e:
    fail("could not reach Keycloak at %s -- is the container running?" % KC_URL, e)
token = tok["access_token"]

# Resolve the group once. A typo here would otherwise produce accounts
# that log in and are then denied every application, which looks like a
# ShinyProxy problem rather than a missing group.
try:
    _, _, groups = call("GET", "/admin/realms/%s/groups?search=%s"
                        % (REALM, urllib.parse.quote(GROUP)), token=token)
except urllib.error.HTTPError as e:
    fail("could not list groups in realm '%s'" % REALM, e)
match = [g for g in (groups or []) if g["name"] == GROUP]
if not match:
    fail("realm '%s' has no group '%s' -- was omicsapp-realm.json imported?"
         % (REALM, GROUP))
group_id = match[0]["id"]

# No look-alike characters: these are read aloud and typed by hand.
ALPHABET = "".join(c for c in string.ascii_letters + string.digits
                   if c not in "O0oIl1")

created = []
for email in emails:
    try:
        status, headers, _ = call("POST", "/admin/realms/%s/users" % REALM, token=token,
                                  body={"username": email, "email": email,
                                        "enabled": True, "emailVerified": False})
    except urllib.error.HTTPError as e:
        if e.code == 409:
            print("skip   %s -- already exists" % email)
            continue
        fail("could not create %s" % email, e)

    # Keycloak returns the new account's location; the last segment is
    # the sub, which is what everything downstream is keyed on.
    location = headers.get("Location", "")
    sub = location.rstrip("/").rsplit("/", 1)[-1]
    if len(sub) != 36:
        fail("Keycloak did not return a usable id for %s (got %r)" % (email, sub))

    # Before the group, so that no account can start a container until
    # the directory Docker would otherwise create as root exists.
    made = []
    roots = [(USERS_ROOT, "projects")]
    if RAW_ROOT:
        roots.append((RAW_ROOT, "archived uploads"))
    for root, label in roots:
        path = os.path.join(root, sub)
        if os.path.exists(path):
            print("note   %s already exists; left untouched" % path)
            continue
        os.makedirs(path, mode=0o700)
        os.chown(path, UID, GID)
        os.chmod(path, 0o700)   # makedirs applies the umask; chmod does not
        made.append("%s (%s)" % (path, label))

    pw = "".join(secrets.choice(ALPHABET) for _ in range(14))
    try:
        # temporary=True is what makes Keycloak demand a new password at
        # first login. It is the whole reason the admin never needs to
        # know anyone's real password.
        call("PUT", "/admin/realms/%s/users/%s/reset-password" % (REALM, sub), token=token,
             body={"type": "password", "value": pw, "temporary": True})
        call("PUT", "/admin/realms/%s/users/%s/groups/%s" % (REALM, sub, group_id), token=token)
    except urllib.error.HTTPError as e:
        fail("created %s but could not finish setting it up" % email, e)

    created.append((email, pw, sub))
    print("create %s" % email)
    for m in made:
        print("       %s" % m)

if not created:
    sys.exit(0)

fresh = not os.path.exists(SECRETS)
with open(os.open(SECRETS, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600), "a") as fh:
    if fresh:
        fh.write("%-34s  %-14s  %s\n" % ("ACCOUNT", "TEMP PASSWORD", "STORAGE ID"))
    for email, pw, sub in created:
        fh.write("%-34s  %-14s  %s\n" % (email, pw, sub))

print()
print("%d account(s) created. Temporary passwords: %s" % (len(created), SECRETS))
PY

# Under sudo the handout file belongs to the person who ran it.
if [ -n "${SUDO_USER:-}" ] && [ -f "$SECRETS" ]; then
    chown "$SUDO_USER" "$SECRETS" 2>/dev/null || true
fi

unset KC_ADMIN_PASSWORD

cat <<DONE

Hand each person their address and temporary password. Keycloak will
require a new one at first login, so what is in that file stops working
as soon as they have chosen their own -- which is why you never need to
know it.

Once they are distributed:
    shred -u ${SECRETS}

Nothing here needs a ShinyProxy restart. Accounts live in Keycloak, so
nobody's running analysis is interrupted by adding a colleague.
DONE
