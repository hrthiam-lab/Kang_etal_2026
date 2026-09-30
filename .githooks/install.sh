#!/bin/sh
set -eu
guard_root=$(git rev-parse --show-toplevel)
cd "$guard_root"
guard_existing=$(git config --get core.hooksPath || true)
if [ -n "$guard_existing" ] && [ "$guard_existing" != '.githooks' ]; then
    printf '%s\n' "Existing hooksPath $guard_existing must be reviewed before installation." >&2
    exit 1
fi
if [ -z "$guard_existing" ]; then
    for guard_hook_name in commit-msg pre-push; do
        guard_old_hook=$(git rev-parse --git-path "hooks/$guard_hook_name")
        if [ -f "$guard_old_hook" ]; then
            printf '%s\n' "Existing hook $guard_old_hook must be integrated, not overwritten." >&2
            exit 1
        fi
    done
fi
chmod +x .githooks/commit-msg .githooks/pre-push .githooks/*.sh
git config --local user.name minwookang-stanford
git config --local user.email mwkang@stanford.edu
git config --local user.useConfigOnly true
git config --local pull.ff only
git config --local core.hooksPath .githooks
sh .githooks/check-history.sh
printf '%s\n' 'Contributor guards installed for this clone. Repeat installation in each new clone.'
