# Repository instructions

- Write code comments in English.
- Use `minwookang-stanford <mwkang@stanford.edu>` for Minwoo's commits. Do not use his personal GitHub account or personal email for attribution.
- Never add AI agents, models, or vendors as authors/coauthors, or append AI attribution to commit messages or pull request descriptions.
- Install the repository's Git guards with `sh .githooks/install.sh` in each new clone. Do not disable hooks or use `--no-verify` in this repository.
- Run `sh .githooks/check-history.sh HEAD` before pushing. New authors/coauthors require the owner's approval and an intentional update to `.githooks/allowed-authors`.
- Never merge or push `backup/*`, `refs/original/*`, or any pre-cleanup copy of the history. Do not use `git push --all` or `git push --mirror`.
- After any upstream history rewrite, stop if local and remote histories diverge. Preserve uncommitted work and ask before resynchronizing; do not merge the old history into the cleaned branch, and do not reset or delete user work.
- Do not weaken the guards or modify attribution settings without the owner's explicit approval.
- Tests may deliberately construct invalid commits with hooks disabled only in disposable local fixture repositories. Never push those fixtures to GitHub.

See `CONTRIBUTOR_GUARDS.md` for scope, limitations, and recovery guidance.
