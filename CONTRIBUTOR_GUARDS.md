# Contributor attribution safeguards

## Install in every clone

From Git Bash or another POSIX shell in this repository:

```sh
sh .githooks/install.sh
```

This configures only this clone: Stanford name/email, explicit identity configuration, fast-forward-only pulls, and `core.hooksPath=.githooks`. The installer refuses to replace unrelated hooks. Once installed, normal Git commits and pushes (including clients that invoke Git, such as GitHub Desktop) use the guards.

## What is blocked

- Commit authors and coauthors outside `.githooks/allowed-authors` (currently only `mwkang@stanford.edu`). The GitHub web committer `noreply@github.com` is allowed as a committer, not as an author.
- AI coauthor trailers, even if an approved email is used.
- Pushes whose complete reachable history contains invalid attribution, including old commits resurrected by merges and commits referenced by tags.
- Publishing `backup/*` or `refs/original/*` recovery references. Existing local recovery history is preserved and is not part of the normal `main` check.

The `commit-msg` hook also checks the pending merge parents when Git exposes `MERGE_HEAD`. The pre-push full-history check is the final local barrier, including for commits created through paths that skip commit hooks.

## Check and recover safely

```sh
sh .githooks/check-history.sh HEAD
```

If a guard fails, do not bypass it. Preserve the working tree and local commits, then inspect the reported commit. After an upstream rewrite, use a fresh clone or a reviewed resynchronization procedure that preserves local work. An ordinary merge with a pre-cleanup clone can restore the removed authors. Do not blindly reset, merge, or force-push.

To add a genuine collaborator, ask the owner to approve their verified commit email before updating `allowed-authors`. Keep the special prohibition on AI coauthors.

## Claude and other coding agents

`.claude/settings.json` disables commit and PR attribution; `CLAUDE.md` and `AGENTS.md` prohibit AI author/coauthor entries. Higher-priority local/managed Claude settings can override project settings, so the Git guards also validate the actual commit metadata.

## Limits

Git hooks are local and are not automatically installed by cloning. Install them on every PC and every independent checkout. Web edits, other unprotected clones, and explicit hook bypasses can evade local guards. These changes do not impose branch protection or change collaborators' GitHub permissions. GitHub's cached contributor display can still lag independently of a clean history; a local hook cannot control that server cache.

## Tests

```sh
sh .githooks/test.sh
```

Tests use disposable repositories and a local bare remote, never GitHub. Fixture locations are printed and retained for inspection.
