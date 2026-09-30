#!/bin/sh
set -eu
guard_test_hooks=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
guard_test_dir=$(mktemp -d "${TMPDIR:-/tmp}/kang-contributor-guard.XXXXXX")
printf 'Local test fixtures: %s\n' "$guard_test_dir"
git init -q --bare "$guard_test_dir/remote.git"
git init -q --initial-branch=main "$guard_test_dir/work"
cd "$guard_test_dir/work"
git config user.name minwookang-stanford
git config user.email mwkang@stanford.edu
git config core.hooksPath "$guard_test_hooks"
git remote add test-remote "$guard_test_dir/remote.git"

guard_test_fail() {
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

guard_expect_rejection() {
    guard_test_label=$1
    shift
    if "$@" >"$guard_test_dir/rejection.log" 2>&1; then
        guard_test_fail "$guard_test_label unexpectedly succeeded"
    fi
    if ! grep -q 'Contributor guard:' "$guard_test_dir/rejection.log"; then
        cat "$guard_test_dir/rejection.log" >&2
        guard_test_fail "$guard_test_label failed for an unrelated reason"
    fi
    printf 'PASS: %s blocked\n' "$guard_test_label"
}

git commit -q --allow-empty -m 'Approved Stanford identity'
guard_test_clean=$(git rev-parse HEAD)
printf '%s\n' 'PASS: clean commit allowed'
guard_expect_rejection 'AI coauthor' git commit --allow-empty -m 'Invalid AI attribution' -m 'Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>'
guard_expect_rejection 'Personal author override' git commit --allow-empty --author='Minwoo Kang <minwookang91@gmail.com>' -m 'Invalid author'
guard_expect_rejection 'Personal committer override' env GIT_COMMITTER_EMAIL=minwookang91@gmail.com git commit --allow-empty -m 'Invalid committer'
guard_expect_rejection 'Unapproved human coauthor' git commit --allow-empty -m 'Unapproved collaborator' -m 'Co-Authored-By: Example Human <unapproved@example.com>'
[ "$(git rev-parse HEAD)" = "$guard_test_clean" ] || guard_test_fail 'Rejected commits changed HEAD'

# Deliberately create invalid history only inside this isolated local fixture.
git switch -q -c poisoned
env GIT_AUTHOR_EMAIL=minwookang91@gmail.com GIT_COMMITTER_EMAIL=minwookang91@gmail.com git -c core.hooksPath="$guard_test_dir/no-hooks" commit -q --allow-empty -m 'Historical personal identity'
git -c core.hooksPath="$guard_test_dir/no-hooks" commit -q --allow-empty -m 'Historical AI coauthor' -m 'co-authored-by: Claude <noreply@anthropic.com>'
git -c core.hooksPath="$guard_test_dir/no-hooks" commit -q --allow-empty -m 'Clean tip above invalid history'
guard_test_poisoned=$(git rev-parse HEAD)
git switch -q main
git commit -q --allow-empty -m 'Clean divergent local work'
guard_test_before_merge=$(git rev-parse HEAD)
guard_expect_rejection 'Resurrected history merge' git merge --no-edit poisoned
[ "$(git rev-parse HEAD)" = "$guard_test_before_merge" ] || guard_test_fail 'Rejected merge changed HEAD'
git merge --abort
guard_expect_rejection 'Full history despite clean tip' git push test-remote poisoned:refs/heads/poisoned
git tag old-history "$guard_test_poisoned"
guard_expect_rejection 'Tag referencing invalid history' git push test-remote refs/tags/old-history
git branch backup/recovered main
guard_expect_rejection 'Backup publication' git push test-remote backup/recovered:refs/heads/backup/recovered

git push -q test-remote main
printf '%s\n' 'PASS: clean push allowed'
[ "$(git --git-dir="$guard_test_dir/remote.git" for-each-ref --format='%(refname)' | wc -l | tr -d ' ')" = 1 ] || guard_test_fail 'Rejected references reached local remote'
sh "$guard_test_hooks/check-history.sh" main
printf '%s\n' 'All contributor guard integration tests passed.'
