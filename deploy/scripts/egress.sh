#!/usr/bin/env bash
#
# Keep the app containers off the network: they may answer connections,
# never open them.
#
#   sudo deploy/scripts/egress.sh apply     # add the rules (idempotent)
#   sudo deploy/scripts/egress.sh check     # exit 0 only if they are in place
#   sudo deploy/scripts/egress.sh remove    # take them out again
#
# The app processes user-supplied data with 224 R packages, and nothing
# it does needs a connection out: gene sets are in the image or on the
# pre-warmed volume, and ShinyProxy is the one that connects *in*. A
# container that can reach the internet, the LAN or the host's services
# is a container a compromised package can exfiltrate data from, or use
# as a foothold, so every new connection from the app network is dropped:
#
#   DOCKER-USER  -i <bridge> NEW -> DROP  forwarded traffic: the internet,
#                                          the LAN, other Docker networks
#                                          (Keycloak and its database
#                                          included)
#   INPUT        -i <bridge> NEW -> DROP  the host's own services (sshd,
#                                          nginx, a colleague's :8080)
#
# Replies on connections ShinyProxy opened are ESTABLISHED, not NEW, so
# the app keeps working. Not `docker network create --internal`: that
# also removes published ports, which is how ShinyProxy reaches a
# container when it runs on the host (internal-networking: false).
#
# The bridge is found from the network's fixed name option (README,
# step 5: -o com.docker.network.bridge.name=br-omicsapp), so the rules
# survive the network being recreated. Docker never flushes DOCKER-USER;
# a reboot does, which is what deploy/systemd/omicsapp-egress.service is
# for.
#
# Environment: NETWORK (default omicsapp-net), DOCKER, IPTABLES,
# IP6TABLES (for testing).

set -euo pipefail

NETWORK="${NETWORK:-omicsapp-net}"
DOCKER="${DOCKER:-docker}"
IPTABLES="${IPTABLES:-iptables}"
IP6TABLES="${IP6TABLES:-ip6tables}"
TAG="omicsapp-egress"

die() { echo "egress.sh: $*" >&2; exit 1; }

action="${1:-}"
case "$action" in
    apply|check|remove) ;;
    *) sed -n '2,8p' "$0"; exit 2 ;;
esac

bridge="$("$DOCKER" network inspect -f '{{index .Options "com.docker.network.bridge.name"}}' "$NETWORK" 2>/dev/null || true)"
if [ -z "$bridge" ] || [ "$bridge" = "<no value>" ]; then
    id="$("$DOCKER" network inspect -f '{{.Id}}' "$NETWORK" 2>/dev/null || true)"
    [ -n "$id" ] || die "no Docker network '$NETWORK'. Create it as in deploy/README.md, step 5."
    # Docker's default name for a user-defined bridge.
    bridge="br-${id:0:12}"
    echo "egress.sh: '$NETWORK' has no fixed bridge name; using $bridge, which changes if the network is recreated." >&2
fi

# "<chain> <match...>", one rule per line.
rules() {
    printf '%s\n' \
        "DOCKER-USER -i $bridge -m conntrack --ctstate NEW,INVALID -m comment --comment $TAG -j DROP" \
        "INPUT -i $bridge -m conntrack --ctstate NEW,INVALID -m comment --comment $TAG -j DROP"
}

missing=0
for ipt in "$IPTABLES" "$IP6TABLES"; do
    command -v "$ipt" >/dev/null 2>&1 || { [ "$ipt" = "$IPTABLES" ] && die "$ipt not found"; continue; }
    if ! "$ipt" -n -L DOCKER-USER >/dev/null 2>&1; then
        [ "$ipt" = "$IP6TABLES" ] && continue   # IPv6 not enabled in Docker: nothing to do
        die "no DOCKER-USER chain in $ipt. Is Docker running with its iptables backend? (With firewall-backend=nftables, add the equivalent rule to Docker's nftables setup instead.)"
    fi
    while IFS= read -r rule; do
        chain="${rule%% *}"
        spec="${rule#* }"
        # shellcheck disable=SC2086
        if "$ipt" -C "$chain" $spec 2>/dev/null; then present=1; else present=0; fi
        case "$action" in
            apply)
                # shellcheck disable=SC2086
                [ "$present" = 1 ] || "$ipt" -I "$chain" 1 $spec
                ;;
            remove)
                # shellcheck disable=SC2086
                while "$ipt" -C "$chain" $spec 2>/dev/null; do "$ipt" -D "$chain" $spec; done
                ;;
            check)
                if [ "$present" = 0 ]; then
                    echo "missing ($ipt): $rule" >&2
                    missing=$((missing + 1))
                fi
                ;;
        esac
    done < <(rules)
done

case "$action" in
    apply)  echo "egress: new connections from $bridge ($NETWORK) are dropped" ;;
    remove) echo "egress: rules for $bridge removed" ;;
    check)
        [ "$missing" = 0 ] || die "$missing rule(s) missing; run: sudo $0 apply"
        echo "egress: in place for $bridge ($NETWORK)"
        ;;
esac
