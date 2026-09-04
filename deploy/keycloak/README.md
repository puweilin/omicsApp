# Keycloak

Keycloak holds the accounts. ShinyProxy no longer stores a password
anywhere: it redirects to Keycloak, Keycloak decides, and ShinyProxy
trusts the result.

That indirection buys one thing — **users can change their own password
and their own username** without an admin editing a file and restarting
a service that would drop everyone else's running analysis.

## The one design decision to understand

Storage is keyed on the Keycloak `sub`, not on a name:

```
/srv/omicsapp/users/f47ac10b-58cc-4372-a567-0e02b2c3d479/
```

`sub` is a UUID fixed when the account is created and never changes.
Usernames and email addresses can change, and the day someone renames
themselves a name-keyed directory would stop being theirs — the app
would open an empty workspace and their projects would still be sitting
in a directory nothing points at.

So the three settings below are a set. Changing one without the others
loses data:

| where | setting | value |
|---|---|---|
| `application.yml` | `proxy.openid.username-attribute` | `sub` |
| `omicsapp-realm.json` | `editUsernameAllowed` | `true` |
| `add_user.sh` | directory name | the `sub` |

The cost is that `ls /srv/omicsapp/users` is unreadable.
`deploy/scripts/list_users.sh` maps it back, and `--link` maintains a
`by-name/` tree of symlinks for when you are chasing a quota.

## First-time setup

Order matters: HTTPS before Keycloak, because Keycloak refuses to serve
a non-local hostname over plain HTTP; Keycloak before ShinyProxy,
because ShinyProxy fails to start if it cannot fetch the signing keys.

**1. Certificate and nginx.** Both commands are in
`deploy/nginx/omicsapp.conf`, including the `update-ca-certificates`
step. Do not skip that one — see *Things that will bite you*.

**2. Secrets.**

```bash
cd deploy/keycloak && cp .env.template .env && chmod 600 .env
```

Generate both values on the server; the template shows how.

**3. Start it.**

```bash
sudo docker compose up -d
sudo docker compose logs -f keycloak      # wait for "Listening on"
```

The realm is imported on this first start, while the database is empty.
After that the database is authoritative and editing
`omicsapp-realm.json` does nothing.

**4. The client secret — the one manual step.**

Admin console → `omicsapp` realm → Clients → `shinyproxy` → Credentials
→ Regenerate. Copy it into `client-secret:` in `application.yml`.

It is deliberately not in this repository. A secret that is never
written to a tracked file cannot be committed by accident, and even in a
private repository a secret committed once stays in the history.

**5. Point ShinyProxy at it.**

```bash
sudo systemctl restart shinyproxy
sudo journalctl -u shinyproxy -n 50
```

**6. Your own account.**

```bash
sudo deploy/scripts/add_user.sh you@example.com
```

Then add yourself to `admins` in the Keycloak console — `add_user.sh`
grants `lab` only, because handing out administrative access should
take a deliberate second action.

## Adding people

```bash
sudo deploy/scripts/add_user.sh alice@example.com bob@example.com
```

Per account this creates the Keycloak user with the email as the
username, generates a **temporary** password, creates the storage
directory named after the `sub` and owned by the container uid, and adds
the account to `lab`.

Hand out the address and the temporary password. Keycloak demands a new
password at first login, so what you handed out stops working as soon as
they have chosen their own — **you never learn anyone's real password.**

No restart. Nobody's analysis is interrupted by a colleague joining.

Delete the handout file once distributed: `shred -u omicsapp-passwords.txt`

## What users can do themselves

At `https://192.168.51.52/auth/realms/omicsapp/account`:

- change their password
- change their username
- change their display name and email
- enrol a second factor

None of it touches the server, and none of it moves their data.

## Realm settings, and why

Stored in `omicsapp-realm.json`. JSON has no comments and Keycloak
rejects unknown keys outright — `FAIL_ON_UNKNOWN_PROPERTIES` is left at
Jackson's default and `RealmRepresentation` carries no
`@JsonIgnoreProperties` — so a `_comment` key would fail the import
rather than be ignored. Hence this table.

| setting | value | why |
|---|---|---|
| `registrationAllowed` | `false` | Accounts come from `add_user.sh`, which also creates the storage directory. A self-registered account would have none, and Docker would create one owned by root that the container cannot write to. |
| `editUsernameAllowed` | `true` | The requirement. Safe only because storage is keyed on `sub`. |
| `loginWithEmailAllowed` | `true` | Accounts start with the email as the username; someone who renames themselves can still sign in the way they are used to. |
| `resetPasswordAllowed` | `false` | Sends mail, and there is no SMTP server. Enabled without one, the forgot-password link fails *after* the user has been told it was sent. Resets are done by an admin. |
| `verifyEmail` | `false` | Same reason. |
| `bruteForceProtected` | `true` | Locks an account for a minute after ten failures. |
| `passwordPolicy` | `length(12)` | Applies to what users choose, not to the generated temporary ones. |
| groups mapper on `shinyproxy` | `groups`, `full.path: false` | **The one that gets missed.** ShinyProxy's `roles-claim` reads it from the ID token and the userinfo response, so all three claim flags are on. Without the mapper every login succeeds and every user is then denied, with nothing in the log connecting the two. |

The redirect URI is not a choice: ShinyProxy's `OpenIDConfiguration.java`
sets `REG_ID = "shinyproxy"` and the template
`{baseUrl}/login/oauth2/code/{registrationId}`.

## Things that will bite you

**Java does not use the system trust store.** ShinyProxy fetches tokens
and signing keys from Keycloak directly — server to server, not through
the browser — and with a self-signed certificate it will refuse. The
browser login looks perfect and the callback then fails with a TLS error
inside a stack trace. `update-ca-certificates` is what fixes it, and the
certificate needs `subjectAltName=IP:192.168.51.52`, because Java
validates against the SAN list and ignores CN.

**`KC_PROXY_HEADERS` missing means 403.** Not a redirect loop, not a
certificate error — a flat 403 whose message says nothing about proxies.

**The bootstrap admin is read once.** `KC_ADMIN_PASSWORD` in `.env` only
applies while the database is empty. Afterwards the account lives in
Postgres; change it in the console.

**Backing up `/srv/omicsapp/users` is not enough.** Without
`/srv/omicsapp/keycloak-db` it is a backup of files nobody can log in to
reach. Both are on the HDD array so one job covers them, but the job has
to name both.

**A deleted user leaves their directory behind.** Nothing removes it
automatically — it holds the only copy of that person's work.
`list_users.sh` reports these as orphans.
